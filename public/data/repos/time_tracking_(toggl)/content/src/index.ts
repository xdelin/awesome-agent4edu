#!/usr/bin/env node
import { Server } from "@modelcontextprotocol/sdk/server/index.js";
import { StdioServerTransport } from "@modelcontextprotocol/sdk/server/stdio.js";
import {
  ListToolsRequestSchema,
  CallToolRequestSchema,
  ErrorCode,
  McpError,
} from "@modelcontextprotocol/sdk/types.js";
import axios, { AxiosInstance } from "axios";
import { z } from "zod";
import dotenv from "dotenv";

dotenv.config();

const TOGGL_API_KEY = process.env.TOGGL_API_KEY;

if (!TOGGL_API_KEY) {
  console.error("Error: TOGGL_API_KEY environment variable is required");
  console.error("Get your API key from: https://track.toggl.com/profile");
  process.exit(1);
}

class TogglClient {
  private api: AxiosInstance;
  private workspaceId?: number;

  constructor(apiKey: string) {
    this.api = axios.create({
      baseURL: "https://api.track.toggl.com/api/v9",
      auth: {
        username: apiKey,
        password: "api_token",
      },
      headers: {
        "Content-Type": "application/json",
      },
    });
  }

  async getWorkspaceId(): Promise<number> {
    if (this.workspaceId) return this.workspaceId;

    const response = await this.api.get("/me");
    this.workspaceId = response.data.default_workspace_id;
    return this.workspaceId!;
  }

  async getCurrentEntry() {
    const response = await this.api.get("/me/time_entries/current");
    return response.data;
  }

  async startTimer(description: string, projectId?: number) {
    const workspaceId = await this.getWorkspaceId();
    const data = {
      created_with: "toggl-mcp",
      description,
      workspace_id: workspaceId,
      project_id: projectId,
      start: new Date().toISOString(),
      duration: -1,
    };

    const response = await this.api.post(`/workspaces/${workspaceId}/time_entries`, data);
    return response.data;
  }

  async stopTimer() {
    const current = await this.getCurrentEntry();
    if (!current) {
      throw new Error("No timer is currently running");
    }

    const workspaceId = await this.getWorkspaceId();
    const response = await this.api.patch(
      `/workspaces/${workspaceId}/time_entries/${current.id}/stop`
    );
    return response.data;
  }

  async getTodayEntries() {
    const workspaceId = await this.getWorkspaceId();
    const today = new Date();
    today.setHours(0, 0, 0, 0);

    const response = await this.api.get("/me/time_entries", {
      params: {
        start_date: today.toISOString(),
        end_date: new Date().toISOString(),
      },
    });

    return response.data;
  }

  async getProjects() {
    const workspaceId = await this.getWorkspaceId();
    const response = await this.api.get(`/workspaces/${workspaceId}/projects`);
    return response.data;
  }

  async deleteEntry(entryId: number) {
    const workspaceId = await this.getWorkspaceId();
    await this.api.delete(`/workspaces/${workspaceId}/time_entries/${entryId}`);
    return { success: true, message: `Deleted entry ${entryId}` };
  }

  async getEntriesForDateRange(startDate: Date, endDate: Date) {
    const response = await this.api.get("/me/time_entries", {
      params: {
        start_date: startDate.toISOString(),
        end_date: endDate.toISOString(),
      },
    });
    return response.data;
  }

  async getWeeklyEntries(weekOffset: number = 0) {
    // Calculate Monday of the target week
    const today = new Date();
    const dayOfWeek = today.getDay();
    const daysToMonday = dayOfWeek === 0 ? -6 : 1 - dayOfWeek;

    const monday = new Date(today);
    monday.setDate(today.getDate() + daysToMonday + (weekOffset * 7));
    monday.setHours(0, 0, 0, 0);

    const sunday = new Date(monday);
    sunday.setDate(monday.getDate() + 6);
    sunday.setHours(23, 59, 59, 999);

    const entries = await this.getEntriesForDateRange(monday, sunday);

    // Calculate summaries
    let totalSeconds = 0;
    const dailyTotals: Record<string, number> = {};
    const projectTotals: Record<string, number> = {};

    for (const entry of entries) {
      let duration: number;
      if (entry.duration < 0) {
        // Timer is running - calculate elapsed time from start
        const startTime = new Date(entry.start).getTime();
        const now = Date.now();
        duration = Math.floor((now - startTime) / 1000);
      } else {
        // Timer is stopped - use the stored duration
        duration = Math.abs(entry.duration || 0);
      }
      totalSeconds += duration;

      // Daily breakdown
      const date = entry.start.split('T')[0];
      dailyTotals[date] = (dailyTotals[date] || 0) + duration;

      // Project breakdown
      const projectName = entry.project_name || "No Project";
      projectTotals[projectName] = (projectTotals[projectName] || 0) + duration;
    }

    return {
      week_starting: monday.toISOString().split('T')[0],
      week_ending: sunday.toISOString().split('T')[0],
      total_hours: Math.round(totalSeconds / 3600 * 100) / 100,
      daily_breakdown: Object.fromEntries(
        Object.entries(dailyTotals).map(([date, seconds]) =>
          [date, Math.round(seconds / 3600 * 100) / 100]
        ).sort()
      ),
      project_breakdown: Object.fromEntries(
        Object.entries(projectTotals)
          .map(([project, seconds]) =>
            [project, Math.round(seconds / 3600 * 100) / 100]
          )
          .sort((a, b) => (b[1] as number) - (a[1] as number))
      ),
      entry_count: entries.length,
      entries: entries
    };
  }
}

const StartTimerSchema = z.object({
  description: z.string().describe("Description of the task"),
  project_name: z.string().optional().describe("Optional project name"),
});

const DeleteEntrySchema = z.object({
  entry_id: z.number().describe("ID of the time entry to delete"),
});

const WeeklyEntriesSchema = z.object({
  week_offset: z.number().optional().default(0).describe("Number of weeks before current week (0 = current, -1 = last week)"),
});

async function main() {
  const client = new TogglClient(TOGGL_API_KEY!);
  const server = new Server({
    name: "toggl-mcp",
    version: "0.1.0",
  }, {
    capabilities: {
      tools: {},
    },
  });

  server.setRequestHandler(ListToolsRequestSchema, async () => ({
    tools: [
      {
        name: "toggl_start",
        description: "Start a new time tracking entry",
        inputSchema: {
          type: "object",
          properties: {
            description: {
              type: "string",
              description: "Description of the task",
            },
            project_name: {
              type: "string",
              description: "Optional project name",
            },
          },
          required: ["description"],
        },
      },
      {
        name: "toggl_stop",
        description: "Stop the currently running timer",
        inputSchema: {
          type: "object",
          properties: {},
        },
      },
      {
        name: "toggl_current",
        description: "Get the currently running time entry",
        inputSchema: {
          type: "object",
          properties: {},
        },
      },
      {
        name: "toggl_today",
        description: "Get today's time entries with total duration",
        inputSchema: {
          type: "object",
          properties: {},
        },
      },
      {
        name: "toggl_projects",
        description: "List all projects in the workspace",
        inputSchema: {
          type: "object",
          properties: {},
        },
      },
      {
        name: "toggl_delete",
        description: "Delete a time entry by ID",
        inputSchema: {
          type: "object",
          properties: {
            entry_id: {
              type: "number",
              description: "ID of the time entry to delete",
            },
          },
          required: ["entry_id"],
        },
      },
      {
        name: "toggl_weekly",
        description: "Get weekly time tracking summary with total hours, daily/project breakdowns",
        inputSchema: {
          type: "object",
          properties: {
            week_offset: {
              type: "number",
              description: "Number of weeks before current week (0 = current, -1 = last week)",
            },
          },
        },
      },
      {
        name: "toggl_last_week",
        description: "Get last week's time tracking summary",
        inputSchema: {
          type: "object",
          properties: {},
        },
      },
    ],
  }));

  server.setRequestHandler(CallToolRequestSchema, async (request: any) => {
    try {
      const { name, arguments: args } = request.params;

      switch (name) {
        case "toggl_start": {
          const { description, project_name } = StartTimerSchema.parse(args);

          let projectId: number | undefined;
          if (project_name) {
            const projects = await client.getProjects();
            const project = projects.find(
              (p: any) => p.name.toLowerCase() === project_name.toLowerCase()
            );
            projectId = project?.id;
          }

          const entry = await client.startTimer(description, projectId);
          return {
            content: [
              {
                type: "text",
                text: `Started timer: "${description}"${
                  project_name ? ` on project "${project_name}"` : ""
                }\nEntry ID: ${entry.id}`,
              },
            ],
          };
        }

        case "toggl_stop": {
          const entry = await client.stopTimer();
          const duration = Math.abs(entry.duration);
          const hours = Math.floor(duration / 3600);
          const minutes = Math.floor((duration % 3600) / 60);

          return {
            content: [
              {
                type: "text",
                text: `Stopped timer: "${entry.description}"\nDuration: ${hours}h ${minutes}m`,
              },
            ],
          };
        }

        case "toggl_current": {
          const current = await client.getCurrentEntry();

          if (!current) {
            return {
              content: [
                {
                  type: "text",
                  text: "No timer is currently running",
                },
              ],
            };
          }

          // For running timers, calculate duration from start time
          let duration: number;
          if (current.duration < 0) {
            // Timer is running - calculate elapsed time from start
            const startTime = new Date(current.start).getTime();
            const now = Date.now();
            duration = Math.floor((now - startTime) / 1000);
          } else {
            // Timer is stopped - use the stored duration
            duration = Math.abs(current.duration);
          }

          const hours = Math.floor(duration / 3600);
          const minutes = Math.floor((duration % 3600) / 60);

          return {
            content: [
              {
                type: "text",
                text: `Currently tracking: "${current.description}"\nRunning for: ${hours}h ${minutes}m`,
              },
            ],
          };
        }

        case "toggl_today": {
          const entries = await client.getTodayEntries();

          if (!entries || entries.length === 0) {
            return {
              content: [
                {
                  type: "text",
                  text: "No time entries today",
                },
              ],
            };
          }

          const totalSeconds = entries.reduce((sum: number, entry: any) => {
            let duration: number;
            if (entry.duration < 0) {
              // Timer is running - calculate elapsed time from start
              const startTime = new Date(entry.start).getTime();
              const now = Date.now();
              duration = Math.floor((now - startTime) / 1000);
            } else {
              // Timer is stopped - use the stored duration
              duration = Math.abs(entry.duration);
            }
            return sum + duration;
          }, 0);

          const totalHours = Math.floor(totalSeconds / 3600);
          const totalMinutes = Math.floor((totalSeconds % 3600) / 60);

          const entryList = entries
            .map((e: any) => {
              let duration: number;
              if (e.duration < 0) {
                // Timer is running - calculate elapsed time from start
                const startTime = new Date(e.start).getTime();
                const now = Date.now();
                duration = Math.floor((now - startTime) / 1000);
              } else {
                // Timer is stopped - use the stored duration
                duration = Math.abs(e.duration);
              }
              const h = Math.floor(duration / 3600);
              const m = Math.floor((duration % 3600) / 60);
              const runningIndicator = e.duration < 0 ? " (running)" : "";
              return `- ${e.description || "No description"} (${h}h ${m}m)${runningIndicator}`;
            })
            .join("\n");

          return {
            content: [
              {
                type: "text",
                text: `Today's entries:\n${entryList}\n\nTotal: ${totalHours}h ${totalMinutes}m`,
              },
            ],
          };
        }

        case "toggl_projects": {
          const projects = await client.getProjects();

          if (!projects || projects.length === 0) {
            return {
              content: [
                {
                  type: "text",
                  text: "No projects found",
                },
              ],
            };
          }

          const projectList = projects
            .map((p: any) => `- ${p.name} (ID: ${p.id})`)
            .join("\n");

          return {
            content: [
              {
                type: "text",
                text: `Projects:\n${projectList}`,
              },
            ],
          };
        }

        case "toggl_delete": {
          const { entry_id } = DeleteEntrySchema.parse(args);
          const result = await client.deleteEntry(entry_id);

          return {
            content: [
              {
                type: "text",
                text: result.message,
              },
            ],
          };
        }

        case "toggl_weekly": {
          const { week_offset } = WeeklyEntriesSchema.parse(args);
          const summary = await client.getWeeklyEntries(week_offset);

          const weekLabel = week_offset === 0 ? "Current" : week_offset === -1 ? "Last" : `${Math.abs(week_offset)} weeks ago`;

          let text = `${weekLabel} Week Summary (${summary.week_starting} to ${summary.week_ending})\n`;
          text += `Total: ${summary.total_hours} hours\n\n`;

          if (Object.keys(summary.daily_breakdown).length > 0) {
            text += "Daily Breakdown:\n";
            for (const [date, hours] of Object.entries(summary.daily_breakdown)) {
              text += `  ${date}: ${hours}h\n`;

              // Show entries for this day
              const dayEntries = summary.entries.filter((e: any) =>
                e.start && e.start.startsWith(date)
              );

              if (dayEntries.length > 0) {
                const entrySummary = dayEntries
                  .map((e: any) => {
                    let duration: number;
                    if (e.duration < 0) {
                      // Timer is running - calculate elapsed time from start
                      const startTime = new Date(e.start).getTime();
                      const now = Date.now();
                      duration = Math.floor((now - startTime) / 1000);
                    } else {
                      // Timer is stopped - use the stored duration
                      duration = Math.abs(e.duration || 0);
                    }
                    const h = Math.floor(duration / 3600);
                    const m = Math.floor((duration % 3600) / 60);
                    const timeStr = h > 0 ? `${h}h${m > 0 ? ` ${m}m` : ''}` : `${m}m`;
                    const runningIndicator = e.duration < 0 ? " (running)" : "";
                    return `    • ${e.description || "No description"} (${timeStr})${runningIndicator}`;
                  })
                  .join("\n");
                text += entrySummary + "\n";
              }
            }
            text += "\n";
          }

          if (Object.keys(summary.project_breakdown).length > 0) {
            text += "Project Breakdown:\n";
            for (const [project, hours] of Object.entries(summary.project_breakdown)) {
              text += `  ${project}: ${hours}h\n`;
            }
          }

          text += `\nTotal entries: ${summary.entry_count}`;

          return {
            content: [
              {
                type: "text",
                text,
              },
            ],
          };
        }

        case "toggl_last_week": {
          const summary = await client.getWeeklyEntries(-1);

          let text = `Last Week Summary (${summary.week_starting} to ${summary.week_ending})\n`;
          text += `Total: ${summary.total_hours} hours\n\n`;

          if (Object.keys(summary.daily_breakdown).length > 0) {
            text += "Daily Breakdown:\n";
            for (const [date, hours] of Object.entries(summary.daily_breakdown)) {
              text += `  ${date}: ${hours}h\n`;

              // Show entries for this day
              const dayEntries = summary.entries.filter((e: any) =>
                e.start && e.start.startsWith(date)
              );

              if (dayEntries.length > 0) {
                const entrySummary = dayEntries
                  .map((e: any) => {
                    let duration: number;
                    if (e.duration < 0) {
                      // Timer is running - calculate elapsed time from start
                      const startTime = new Date(e.start).getTime();
                      const now = Date.now();
                      duration = Math.floor((now - startTime) / 1000);
                    } else {
                      // Timer is stopped - use the stored duration
                      duration = Math.abs(e.duration || 0);
                    }
                    const h = Math.floor(duration / 3600);
                    const m = Math.floor((duration % 3600) / 60);
                    const timeStr = h > 0 ? `${h}h${m > 0 ? ` ${m}m` : ''}` : `${m}m`;
                    const runningIndicator = e.duration < 0 ? " (running)" : "";
                    return `    • ${e.description || "No description"} (${timeStr})${runningIndicator}`;
                  })
                  .join("\n");
                text += entrySummary + "\n";
              }
            }
            text += "\n";
          }

          if (Object.keys(summary.project_breakdown).length > 0) {
            text += "Project Breakdown:\n";
            for (const [project, hours] of Object.entries(summary.project_breakdown)) {
              text += `  ${project}: ${hours}h\n`;
            }
          }

          text += `\nTotal entries: ${summary.entry_count}`;

          return {
            content: [
              {
                type: "text",
                text,
              },
            ],
          };
        }

        default:
          throw new McpError(
            ErrorCode.MethodNotFound,
            `Unknown tool: ${name}`
          );
      }
    } catch (error: any) {
      if (error instanceof z.ZodError) {
        throw new McpError(
          ErrorCode.InvalidParams,
          `Invalid parameters: ${error.errors.map((e) => e.message).join(", ")}`
        );
      }

      if (error.response?.status === 401) {
        throw new McpError(
          ErrorCode.InvalidRequest,
          "Invalid Toggl API key. Check your TOGGL_API_KEY environment variable."
        );
      }

      throw new McpError(
        ErrorCode.InternalError,
        error.message || "An unexpected error occurred"
      );
    }
  });

  const transport = new StdioServerTransport();
  await server.connect(transport);

  console.error("Toggl MCP server running");
}

main().catch((error) => {
  console.error("Fatal error:", error);
  process.exit(1);
});