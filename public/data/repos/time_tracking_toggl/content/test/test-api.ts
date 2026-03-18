#!/usr/bin/env tsx
import axios, { AxiosInstance } from "axios";
import dotenv from "dotenv";

dotenv.config();

const TOGGL_API_KEY = process.env.TOGGL_API_KEY || "d0b1a0cd530218adb643580eb4816e7c";

class TogglTester {
  private api: AxiosInstance;

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

  async testAll() {
    console.log("🧪 Testing Toggl API with your credentials...\n");

    try {
      // Test 1: Get user info and workspace
      console.log("1️⃣ Testing authentication and getting workspace...");
      const meResponse = await this.api.get("/me");
      const workspaceId = meResponse.data.default_workspace_id;
      console.log(`✅ Authenticated as: ${meResponse.data.fullname}`);
      console.log(`   Email: ${meResponse.data.email}`);
      console.log(`   Workspace ID: ${workspaceId}\n`);

      // Test 2: Get projects
      console.log("2️⃣ Testing project listing...");
      const projectsResponse = await this.api.get(`/workspaces/${workspaceId}/projects`);
      const projects = projectsResponse.data || [];
      console.log(`✅ Found ${projects.length} projects:`);
      projects.slice(0, 5).forEach((p: any) => {
        console.log(`   - ${p.name} (ID: ${p.id})`);
      });
      if (projects.length > 5) console.log(`   ... and ${projects.length - 5} more`);
      console.log();

      // Test 3: Get current timer
      console.log("3️⃣ Testing current timer check...");
      const currentResponse = await this.api.get("/me/time_entries/current");
      if (currentResponse.data) {
        console.log(`✅ Current timer running: "${currentResponse.data.description}"`);
      } else {
        console.log("✅ No timer currently running");
      }
      console.log();

      // Test 4: Get today's entries
      console.log("4️⃣ Testing today's entries...");
      const today = new Date();
      today.setHours(0, 0, 0, 0);
      const entriesResponse = await this.api.get("/me/time_entries", {
        params: {
          start_date: today.toISOString(),
          end_date: new Date().toISOString(),
        },
      });
      const entries = entriesResponse.data || [];
      const totalSeconds = entries.reduce((sum: number, entry: any) => {
        return sum + Math.abs(entry.duration || 0);
      }, 0);
      const hours = Math.floor(totalSeconds / 3600);
      const minutes = Math.floor((totalSeconds % 3600) / 60);

      console.log(`✅ Found ${entries.length} entries today (Total: ${hours}h ${minutes}m):`);
      entries.slice(0, 3).forEach((e: any) => {
        const dur = Math.abs(e.duration || 0);
        const h = Math.floor(dur / 3600);
        const m = Math.floor((dur % 3600) / 60);
        console.log(`   - ${e.description || "No description"} (${h}h ${m}m)`);
      });
      if (entries.length > 3) console.log(`   ... and ${entries.length - 3} more`);
      console.log();

      // Test 5: Create and stop a test timer
      console.log("5️⃣ Testing timer creation and stopping...");
      console.log("   Creating test timer...");
      const testEntry = {
        created_with: "toggl-mcp-test",
        description: "MCP Test Timer - Delete Me",
        workspace_id: workspaceId,
        start: new Date().toISOString(),
        duration: -1,
      };

      const createResponse = await this.api.post(
        `/workspaces/${workspaceId}/time_entries`,
        testEntry
      );
      const createdId = createResponse.data.id;
      console.log(`   ✅ Created timer with ID: ${createdId}`);

      // Wait 2 seconds
      await new Promise(resolve => setTimeout(resolve, 2000));

      // Stop it
      console.log("   Stopping test timer...");
      const stopResponse = await this.api.patch(
        `/workspaces/${workspaceId}/time_entries/${createdId}/stop`
      );
      console.log(`   ✅ Stopped timer successfully`);

      // Delete it
      console.log("   Deleting test timer...");
      await this.api.delete(`/workspaces/${workspaceId}/time_entries/${createdId}`);
      console.log(`   ✅ Deleted test timer\n`);

      console.log("🎉 All tests passed! Your Toggl MCP server should work perfectly!");
      console.log("\n📝 Summary:");
      console.log(`   - API Key is valid`);
      console.log(`   - Can access workspace`);
      console.log(`   - Can list projects (${projects.length} found)`);
      console.log(`   - Can check current timer`);
      console.log(`   - Can get today's entries (${entries.length} found)`);
      console.log(`   - Can create/stop/delete timers`);

    } catch (error: any) {
      console.error("❌ Test failed!");
      if (error.response) {
        console.error(`   Status: ${error.response.status}`);
        console.error(`   Message: ${JSON.stringify(error.response.data)}`);
        if (error.response.status === 403) {
          console.error("\n⚠️  Your API key might be invalid or expired.");
          console.error("   Get a new one from: https://track.toggl.com/profile");
        }
      } else {
        console.error(`   Error: ${error.message}`);
      }
      process.exit(1);
    }
  }
}

// Run tests
const tester = new TogglTester(TOGGL_API_KEY);
tester.testAll().catch(console.error);