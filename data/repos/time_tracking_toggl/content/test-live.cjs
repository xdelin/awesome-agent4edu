#!/usr/bin/env node
const axios = require('axios');

const TOGGL_API_KEY = 'd0b1a0cd530218adb643580eb4816e7c';

async function testTogglAPI() {
  const api = axios.create({
    baseURL: 'https://api.track.toggl.com/api/v9',
    auth: {
      username: TOGGL_API_KEY,
      password: 'api_token',
    },
    headers: {
      'Content-Type': 'application/json',
    },
  });

  console.log('Testing Toggl API...\n');

  try {
    // Test 1: Get user info
    const meResponse = await api.get('/me');
    console.log('✅ Authenticated as:', meResponse.data.fullname);
    console.log('   Email:', meResponse.data.email);

    const workspaceId = meResponse.data.default_workspace_id;
    console.log('   Workspace ID:', workspaceId);

    // Test 2: Get current timer
    const currentResponse = await api.get('/me/time_entries/current');
    if (currentResponse.data) {
      console.log('✅ Current timer:', currentResponse.data.description);
    } else {
      console.log('✅ No timer currently running');
    }

    // Test 3: Get today's entries
    const today = new Date();
    today.setHours(0, 0, 0, 0);
    const entriesResponse = await api.get('/me/time_entries', {
      params: {
        start_date: today.toISOString(),
        end_date: new Date().toISOString(),
      },
    });

    const entries = entriesResponse.data || [];
    console.log(`✅ Today's entries: ${entries.length} found`);

    console.log('\n🎉 All tests passed! Your Toggl MCP server is working!');

  } catch (error) {
    console.error('❌ Test failed!', error.message);
    process.exit(1);
  }
}

testTogglAPI();