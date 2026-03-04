#!/usr/bin/env node

/**
 * mcp_server.js - 真正的 MCP Server 入口
 * 用于连接 Claude Desktop App
 */

const { Server } = require("@modelcontextprotocol/sdk/server/index.js");
const { StdioServerTransport } = require("@modelcontextprotocol/sdk/server/stdio.js");
const { CallToolRequestSchema, ListToolsRequestSchema } = require("@modelcontextprotocol/sdk/types.js");
const html2pptx = require('./scripts/html2pptx');
const pptxgen = require('pptxgenjs');
const path = require('path');
const fs = require('fs');

// 1. 初始化 Server
const server = new Server(
  {
    name: "local-pptx-skill",
    version: "1.0.0",
  },
  {
    capabilities: {
      tools: {},
    },
  }
);

// 读取 SKILL.md 作为工具描述
const skillDescription = fs.readFileSync(path.join(__dirname, 'SKILL.md'), 'utf8');

// 2. 定义工具列表 (List Tools)
server.setRequestHandler(ListToolsRequestSchema, async () => {
  return {
    tools: [
      {
        name: "generate_pptx",
        description: skillDescription, // 把 SKILL.md 喂给 Claude
        inputSchema: {
          type: "object",
          properties: {
            html_content: {
              type: "string",
              description: "The HTML content to render into a PowerPoint slide. Must follow the strict CSS rules defined in the description.",
            },
            filename: {
              type: "string",
              description: "The output filename (e.g., 'presentation.pptx').",
            }
          },
          required: ["html_content", "filename"],
        },
      },
    ],
  };
});

// 3. 处理工具调用 (Call Tool)
server.setRequestHandler(CallToolRequestSchema, async (request) => {
  if (request.params.name === "generate_pptx") {
    const { html_content, filename } = request.params.arguments;
    
    try {
      console.error(`[Server] Generating PPT: ${filename}...`); // Log to stderr (so Claude doesn't see it as output)
      
      const pres = new pptxgen();
      pres.layout = 'LAYOUT_16x9';

      // 临时保存 HTML 文件供 html2pptx 读取 (因为原始 html2pptx 依然是读文件的)
      // 注意：这里用了一个临时文件技巧，因为 html2pptx.js 目前设计是读文件的
      const tempHtmlPath = path.join(__dirname, `temp_${Date.now()}.html`);
      fs.writeFileSync(tempHtmlPath, html_content);

      // 调用核心转换逻辑
      await html2pptx(tempHtmlPath, pres);

      // 保存 PPT
      // 使用绝对路径保存到 downloads 或当前目录
      const outputPath = path.resolve(__dirname, filename || 'output.pptx');
      await pres.writeFile(outputPath);

      // 清理临时文件
      fs.unlinkSync(tempHtmlPath);

      return {
        content: [
          {
            type: "text",
            text: `Successfully generated PowerPoint file at: ${outputPath}`,
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: "text",
            text: `Error generating PPT: ${error.message}`,
          },
        ],
        isError: true,
      };
    }
  }
  throw new Error("Tool not found");
});

// 4. 启动通信
const transport = new StdioServerTransport();
server.connect(transport).catch((error) => {
  console.error("Server error:", error);
  process.exit(1);
});
