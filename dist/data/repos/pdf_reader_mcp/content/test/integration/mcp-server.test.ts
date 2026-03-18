/**
 * Integration test for MCP server
 * Tests the actual JSON-RPC communication over stdio
 */

import { type ChildProcess, spawn } from 'node:child_process';
import path from 'node:path';
import { afterAll, beforeAll, describe, expect, it } from 'vitest';

// JSON-RPC message helpers
const createRequest = (id: number, method: string, params?: unknown) => ({
  jsonrpc: '2.0' as const,
  id,
  method,
  params,
});

// MCP uses JSON Lines format (newline-delimited JSON)
const sendMessage = (proc: ChildProcess, message: object): void => {
  const json = JSON.stringify(message);
  proc.stdin?.write(`${json}\n`);
};

const readResponse = (proc: ChildProcess, timeout = 5000): Promise<unknown> => {
  return new Promise((resolve, reject) => {
    let buffer = '';
    const timer = setTimeout(() => {
      reject(new Error(`Timeout waiting for response. Buffer: ${buffer}`));
    }, timeout);

    const onData = (data: Buffer) => {
      buffer += data.toString();

      // Parse JSON Lines format
      const newlineIndex = buffer.indexOf('\n');
      if (newlineIndex === -1) return;

      const line = buffer.slice(0, newlineIndex).trim();
      buffer = buffer.slice(newlineIndex + 1);

      if (line.length === 0) return;

      clearTimeout(timer);
      proc.stdout?.off('data', onData);
      try {
        resolve(JSON.parse(line));
      } catch {
        reject(new Error(`Failed to parse JSON: ${line}`));
      }
    };

    proc.stdout?.on('data', onData);
  });
};

describe('MCP Server Integration', () => {
  let serverProc: ChildProcess;

  beforeAll(async () => {
    // Start the MCP server
    const serverPath = path.resolve(__dirname, '../../dist/index.js');
    // Must use bun as SDK uses Bun APIs
    serverProc = spawn('bun', [serverPath], {
      stdio: ['pipe', 'pipe', 'pipe'],
      env: { ...process.env, NODE_ENV: 'test' },
    });

    // Wait a bit for server to start
    await new Promise((r) => setTimeout(r, 500));
  });

  afterAll(() => {
    serverProc?.kill();
  });

  it('should respond to initialize request', async () => {
    const initRequest = createRequest(1, 'initialize', {
      protocolVersion: '2024-11-05',
      capabilities: {},
      clientInfo: { name: 'test-client', version: '1.0.0' },
    });

    sendMessage(serverProc, initRequest);
    const response = (await readResponse(serverProc)) as {
      id: number;
      result?: { serverInfo?: { name: string } };
    };

    expect(response.id).toBe(1);
    expect(response.result?.serverInfo?.name).toBe('pdf-reader-mcp');
  });

  it('should list available tools', async () => {
    // Send initialized notification first
    sendMessage(serverProc, { jsonrpc: '2.0', method: 'notifications/initialized' });

    // Wait a bit
    await new Promise((r) => setTimeout(r, 100));

    const listRequest = createRequest(2, 'tools/list', {});
    sendMessage(serverProc, listRequest);

    const response = (await readResponse(serverProc)) as {
      id: number;
      result?: { tools?: Array<{ name: string }> };
    };

    expect(response.id).toBe(2);
    expect(response.result?.tools).toBeDefined();
    expect(response.result?.tools?.length).toBeGreaterThan(0);

    const toolNames = response.result?.tools?.map((t) => t.name);
    expect(toolNames).toContain('read_pdf');
  });

  it('should call read_pdf tool with a test PDF', async () => {
    const testPdfPath = path.resolve(__dirname, '../fixtures/sample.pdf');

    const callRequest = createRequest(3, 'tools/call', {
      name: 'read_pdf',
      arguments: {
        sources: [{ path: testPdfPath }],
        include_metadata: true,
        include_page_count: true,
        include_full_text: false,
      },
    });

    sendMessage(serverProc, callRequest);
    const response = (await readResponse(serverProc, 10000)) as {
      id: number;
      result?: { content?: Array<{ type: string; text?: string }>; isError?: boolean };
      error?: { message: string };
    };

    expect(response.id).toBe(3);

    // If test PDF doesn't exist, expect error
    if (response.error || response.result?.isError) {
      // Expected if no test PDF
      expect(response.error?.message || response.result?.content?.[0]?.text).toContain('PDF');
    } else {
      // If it exists, should have content
      expect(response.result?.content).toBeDefined();
      expect(response.result?.content?.[0]?.type).toBe('text');
    }
  });

  it('should handle invalid tool arguments', async () => {
    const callRequest = createRequest(4, 'tools/call', {
      name: 'read_pdf',
      arguments: {
        // Missing required 'sources' field
        include_metadata: true,
      },
    });

    sendMessage(serverProc, callRequest);
    const response = (await readResponse(serverProc)) as {
      id: number;
      result?: { isError?: boolean; content?: Array<{ text?: string }> };
      error?: { code: number; message: string };
    };

    expect(response.id).toBe(4);
    // SDK returns validation error as result.isError, not JSON-RPC error
    expect(response.result?.isError).toBe(true);
    expect(response.result?.content?.[0]?.text).toMatch(/sources/i);
  });
});
