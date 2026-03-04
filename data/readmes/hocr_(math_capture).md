# HOCR MCP Agent
[![Trust Score](https://archestra.ai/mcp-catalog/api/badge/quality/Wooonster/hocr_mcp_server)](https://archestra.ai/mcp-catalog/wooonster__hocr_mcp_server)

### Vue Frontend

Start command:

```bash
cd hocr-vue-client/

npm install
npm run dev
```

Required packages:
```bash
cd hocr-vue-client/

npm install katex
npm install axios
```

### MCP Server

Start command:

```bash
uvicorn mcp_server:app --reload --host 0.0.0.0 --port 8000
````