/**
 * Simple test client to verify MCP server works
 */

import { spawn } from 'node:child_process';
import path from 'node:path';

const serverProc = spawn('bun', ['dist/index.js'], {
  stdio: ['pipe', 'pipe', 'pipe'],
});

// Capture stderr for debugging
serverProc.stderr?.on('data', (data) => {
  console.error('[stderr]', data.toString());
});

const sendMessage = (message: object): void => {
  const json = JSON.stringify(message);
  console.log('\n[send]', json);
  serverProc.stdin?.write(`${json}\n`);
};

serverProc.stdout?.on('data', (data) => {
  const lines = data.toString().split('\n').filter(Boolean);
  for (const line of lines) {
    try {
      const parsed = JSON.parse(line);
      console.log('[recv]', JSON.stringify(parsed, null, 2));
    } catch {
      console.log('[recv raw]', line);
    }
  }
});

// Test sequence
const runTests = async () => {
  // 1. Initialize
  console.log('\n=== Test 1: Initialize ===');
  sendMessage({
    jsonrpc: '2.0',
    id: 1,
    method: 'initialize',
    params: {
      protocolVersion: '2024-11-05',
      capabilities: {},
      clientInfo: { name: 'test-client', version: '1.0.0' },
    },
  });

  await new Promise((r) => setTimeout(r, 500));

  // 2. Send initialized notification
  console.log('\n=== Test 2: Initialized notification ===');
  sendMessage({
    jsonrpc: '2.0',
    method: 'notifications/initialized',
  });

  await new Promise((r) => setTimeout(r, 200));

  // 3. List tools
  console.log('\n=== Test 3: List tools ===');
  sendMessage({
    jsonrpc: '2.0',
    id: 2,
    method: 'tools/list',
    params: {},
  });

  await new Promise((r) => setTimeout(r, 500));

  // 4. Call read_pdf with a test file (expect error since file doesn't exist)
  console.log('\n=== Test 4: Call read_pdf (missing file) ===');
  sendMessage({
    jsonrpc: '2.0',
    id: 3,
    method: 'tools/call',
    params: {
      name: 'read_pdf',
      arguments: {
        sources: [{ path: '/nonexistent/test.pdf' }],
      },
    },
  });

  await new Promise((r) => setTimeout(r, 1000));

  // 5. Call with invalid args
  console.log('\n=== Test 5: Call read_pdf (invalid args) ===');
  sendMessage({
    jsonrpc: '2.0',
    id: 4,
    method: 'tools/call',
    params: {
      name: 'read_pdf',
      arguments: {
        // Missing required 'sources'
        include_metadata: true,
      },
    },
  });

  await new Promise((r) => setTimeout(r, 500));

  // 6. Call with real PDF if exists
  const testPdfPath = path.resolve(process.cwd(), 'test/fixtures/sample.pdf');
  console.log(`\n=== Test 6: Call read_pdf (${testPdfPath}) ===`);
  sendMessage({
    jsonrpc: '2.0',
    id: 5,
    method: 'tools/call',
    params: {
      name: 'read_pdf',
      arguments: {
        sources: [{ path: testPdfPath }],
        include_metadata: true,
        include_page_count: true,
      },
    },
  });

  await new Promise((r) => setTimeout(r, 2000));

  console.log('\n=== All tests complete ===');
  serverProc.kill();
  process.exit(0);
};

runTests();
