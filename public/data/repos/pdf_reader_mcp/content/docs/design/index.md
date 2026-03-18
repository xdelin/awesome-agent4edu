# Design Philosophy

PDF Reader MCP is built on these core principles:

## 1. Performance First

- **Concurrent Processing** - Multiple PDF sources are processed in parallel
- **Efficient Parsing** - Uses pdfjs-dist for reliable, fast PDF parsing
- **Minimal Overhead** - Direct stdio communication with no HTTP overhead
- **Batch Operations** - Process multiple files in a single request

## 2. Comprehensive Extraction

- **Text Extraction** - Full document or specific pages
- **Page Ranges** - Flexible page selection with ranges like "1-5, 10, 15-20"
- **Metadata Access** - Document properties, author, title, dates
- **Image Extraction** - Embedded images as base64-encoded PNG

## 3. Simple Integration

- **Single Tool** - One `read_pdf` tool handles all extraction needs
- **Standard MCP** - Compatible with any MCP client
- **Easy Setup** - One command installation via npx
- **Multiple Clients** - Works with Claude Desktop, Claude Code, Cursor, and more

## 4. Flexible Input

- **Local Files** - Read PDFs from any path on the filesystem
- **Remote URLs** - Download and process PDFs from URLs
- **Mixed Sources** - Combine local and remote files in one request

## 5. Robust Error Handling

- **Graceful Failures** - One failed source doesn't stop others
- **Clear Errors** - Specific error codes and messages
- **Partial Results** - Get results from successful sources even if some fail

## Technical Stack

- **Runtime**: Node.js 22+
- **PDF Parsing**: pdfjs-dist
- **Image Encoding**: pngjs
- **Schema Validation**: Zod
- **MCP SDK**: @sylphx/mcp-server-sdk
- **Build Tool**: bunup
