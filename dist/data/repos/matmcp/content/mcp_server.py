#!/usr/bin/env python3

from typing import Dict, Optional, Any
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
import uvicorn

app = FastAPI(
    title="GNoME MCP Server",
    description="Model Context Protocol server for GNoME materials database",
    version="0.0.1"
)

# Pydantic models for request/response validation
class ClientInfo(BaseModel):
    name: str
    version: str

class InitializeParams(BaseModel):
    protocolVersion: str = Field(..., description="Version of the MCP protocol")
    clientInfo: ClientInfo

class JsonRpcRequest(BaseModel):
    jsonrpc: str = Field(..., pattern="^2\\.0$")
    method: str
    params: Dict[str, Any]
    id: int

class ServerInfo(BaseModel):
    name: str = "GNoME-MCP-Server"
    version: str = "0.0.1"
    capabilities: Dict[str, Any] = Field(
        default_factory=lambda: {
            "materialsDatabase": {
                "name": "GNoME",
                "provider": "Google DeepMind",
                "version": "1.1.0",
                "endpoint": "https://optimade-gnome.odbx.science/v1"
            },
            "supportedQueries": [
                "structure",
                "property",
                "composition"
            ],
            "supportedProperties": [
                "formation_energy",
                "band_gap",
                "decomposition_energy",
                "space_group"
            ]
        }
    )

class InitializeResult(BaseModel):
    serverInfo: ServerInfo
    protocolVersion: str

class JsonRpcResponse(BaseModel):
    jsonrpc: str = "2.0"
    result: Optional[Dict[str, Any]] = None
    error: Optional[Dict[str, Any]] = None
    id: int

@app.post("/mcp")
async def handle_mcp_request(request: JsonRpcRequest) -> JsonRpcResponse:
    """Handle MCP requests."""
    if request.method == "initialize":
        try:
            # Validate initialize parameters
            params = InitializeParams(**request.params)
            
            # Create initialize response
            result = InitializeResult(
                serverInfo=ServerInfo(),
                protocolVersion=params.protocolVersion
            )
            
            return JsonRpcResponse(
                result=result.model_dump(),
                id=request.id
            )
        except Exception as e:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid initialize parameters: {str(e)}"
            )
    else:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported method: {request.method}"
        )

@app.get("/")
async def root():
    """Root endpoint returning basic server information."""
    return {
        "name": "GNoME MCP Server",
        "version": "0.0.1",
        "endpoints": {
            "/mcp": "MCP protocol endpoint",
            "/docs": "API documentation (Swagger UI)",
            "/redoc": "API documentation (ReDoc)"
        }
    }

if __name__ == "__main__":
    uvicorn.run(
        "mcp_server:app",
        host="0.0.0.0",
        port=8000,
        reload=True
    ) 