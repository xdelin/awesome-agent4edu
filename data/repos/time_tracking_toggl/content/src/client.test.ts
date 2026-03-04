import { describe, it, expect, beforeAll } from "vitest";
import axios from "axios";

const TOGGL_API_KEY = process.env.TOGGL_API_KEY || "d0b1a0cd530218adb643580eb4816e7c";

describe("Toggl API Integration", () => {
  let api: any;
  let workspaceId: number;

  beforeAll(() => {
    api = axios.create({
      baseURL: "https://api.track.toggl.com/api/v9",
      auth: {
        username: TOGGL_API_KEY,
        password: "api_token",
      },
      headers: {
        "Content-Type": "application/json",
      },
    });
  });

  it("should authenticate and get user info", async () => {
    const response = await api.get("/me");
    expect(response.status).toBe(200);
    expect(response.data).toHaveProperty("email");
    expect(response.data).toHaveProperty("fullname");
    expect(response.data).toHaveProperty("default_workspace_id");

    workspaceId = response.data.default_workspace_id;
    expect(workspaceId).toBeGreaterThan(0);
  });

  it("should list projects", async () => {
    const meResponse = await api.get("/me");
    workspaceId = meResponse.data.default_workspace_id;

    const response = await api.get(`/workspaces/${workspaceId}/projects`);
    expect(response.status).toBe(200);
    expect(Array.isArray(response.data)).toBe(true);
  });

  it("should get current timer (or null)", async () => {
    const response = await api.get("/me/time_entries/current");
    expect(response.status).toBe(200);
    // Can be null if no timer running
    if (response.data) {
      expect(response.data).toHaveProperty("id");
      expect(response.data).toHaveProperty("workspace_id");
    }
  });

  it("should get today's entries", async () => {
    const today = new Date();
    today.setHours(0, 0, 0, 0);

    const response = await api.get("/me/time_entries", {
      params: {
        start_date: today.toISOString(),
        end_date: new Date().toISOString(),
      },
    });

    expect(response.status).toBe(200);
    expect(Array.isArray(response.data)).toBe(true);
  });

  it("should create, stop and delete a timer", async () => {
    const meResponse = await api.get("/me");
    workspaceId = meResponse.data.default_workspace_id;

    // Create timer
    const createData = {
      created_with: "toggl-mcp-test",
      description: "Unit Test Timer - Auto Delete",
      workspace_id: workspaceId,
      start: new Date().toISOString(),
      duration: -1,
    };

    const createResponse = await api.post(
      `/workspaces/${workspaceId}/time_entries`,
      createData
    );
    expect(createResponse.status).toBe(200);
    expect(createResponse.data).toHaveProperty("id");

    const entryId = createResponse.data.id;

    // Stop timer
    const stopResponse = await api.patch(
      `/workspaces/${workspaceId}/time_entries/${entryId}/stop`
    );
    expect(stopResponse.status).toBe(200);

    // Delete timer
    const deleteResponse = await api.delete(
      `/workspaces/${workspaceId}/time_entries/${entryId}`
    );
    expect(deleteResponse.status).toBe(200);
  });
});