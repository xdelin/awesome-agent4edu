# Manim MCP

![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-blue)
![Manim](https://img.shields.io/badge/Manim-Animation-green)
![API](https://img.shields.io/badge/FastAPI-REST-teal)

A Docker-based environment for creating mathematical animations with [Manim](https://www.manim.community/), featuring both a CLI interface and a web API with Model Context Protocol (MCP) support for AI assistants.

## 📑 Overview

This project provides:

1. **Containerized Manim Environment**: Run Manim in an isolated, reproducible Docker environment
2. **Web API**: Create and manage Manim animations via HTTP requests
3. **MCP Integration**: Direct interaction with AI assistants like Claude
4. **File Management**: Upload scripts and download generated animations

## 🚀 Quick Start

### Prerequisites

- [Docker](https://docs.docker.com/get-docker/) and Docker Compose installed on your system

### Installation

#### Option 1: Use the Prebuilt Image (Recommended)

Simply pull the prebuilt image from Docker Hub:

```bash
docker pull wstcpyt/manim-docker-mcp:latest
```

Then run it with docker-compose:

```bash
docker compose up -d
```

#### Option 2: Build Locally

1. Clone the repository:
   ```bash
   git clone https://github.com/YOUR_USERNAME/manim-docker-mcp.git
   cd manim-docker-mcp
   ```

2. Build the Docker images:
   ```bash
   docker compose build
   ```

### Usage

#### CLI Mode

Create a Python file in the `animations` directory (see example below), then run:

```bash
docker compose run manim -pql animations/example.py ExampleScene
```

#### API Mode

Start the API server:

```bash
docker compose up -d manim-api
```

Access the API documentation at [http://localhost:8000/docs](http://localhost:8000/docs)

## 🎬 Creating Animations

### Basic Example

Create a file `animations/example.py`:

```python
from manim import *

class CircleToSquare(Scene):
    def construct(self):
        circle = Circle()
        circle.set_fill(BLUE, opacity=0.5)

        square = Square()
        square.set_fill(RED, opacity=0.5)

        self.play(Create(circle))
        self.wait()
        self.play(Transform(circle, square))
        self.wait()
```

### Running the Animation

```bash
# CLI mode with preview (-p), low quality (-ql)
docker compose run manim -pql animations/example.py CircleToSquare

# API mode
curl -X POST "http://localhost:8000/run-manim?filepath=/manim/temp/circle_example.py&scene_name=CircleToSquare&quality=low_quality"
```

## 📂 Project Structure

```
manim-docker-mcp/
├── animations/           # Manim animation scripts
├── app/                  # FastAPI application
├── media/                # Generated animations (CLI mode)
├── output/               # Generated animations (API mode)
├── temp/                 # Temporary files
├── uploads/              # Uploaded animation scripts
├── Dockerfile            # Docker image definition
├── docker-compose.yml    # Docker Compose configuration
└── README.md             # This file
```

## 🔧 Configuration

### Quality Settings

| Flag | Resolution | Frame Rate | Best For |
|------|------------|------------|----------|
| `-ql` | 480p | 15fps | Quick previews |
| `-qm` | 720p | 30fps | General use |
| `-qh` | 1080p | 60fps | Presentations |
| `-qk` | 1440p | 60fps | Production videos |

### Other Useful Flags

- `-p`: Preview the output file
- `-t`: Transparent background
- `--save_last_frame`: Render only the last frame
- `-c COLOR`: Set background color

## 🌐 API Documentation

### Core Endpoints

#### List Files
```http
GET /list-files?directory=/manim
```

#### Write File
```http
POST /write-file?filepath=/manim/temp/example.py
```

#### Run Animation
```http
POST /run-manim?filepath=/manim/temp/example.py&scene_name=CircleToSquare
```

#### Download Animation
```http
GET /download-file?filepath=/media/videos/example/480p15/CircleToSquare.mp4
```

Full API documentation is available at the `/docs` endpoint.

## 🤖 AI Assistant Integration (MCP)

This project supports the [Model Context Protocol (MCP)](https://github.com/tadata-org/fastapi_mcp), enabling AI assistants to:

1. Create Manim scripts based on natural language descriptions
2. Run animations and provide download links
3. Browse and manage generated media files

Example MCP session:

```
User: Create an animation showing a circle morphing into a square
AI: I'll create that for you...
```

## 🔍 Advanced Usage

### Custom LaTeX

The container includes a minimal LaTeX installation. Custom LaTeX can be used in animations:

```python
formula = MathTex(r"\int_{a}^{b} f(x) \, dx = F(b) - F(a)")
self.play(Write(formula))
```

### Mounting Custom Directories

Modify the `docker-compose.yml` file to mount additional directories:

```yaml
volumes:
  - ./my_custom_dir:/manim/custom
```

## 🛠️ Troubleshooting

### Common Issues

- **Docker not running**: Make sure Docker daemon is running
- **Permission errors**: The container needs write access to mounted volumes
- **Missing media**: Check the correct output directory (media/ for CLI, output/ for API)

### Getting Help

If you encounter issues:
1. Check the [Manim documentation](https://docs.manim.community/)
2. Search existing GitHub issues
3. Create a new issue with details about your problem

## 📜 License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgements

- [Manim Community](https://www.manim.community/) for the amazing animation engine
- [FastAPI](https://fastapi.tiangolo.com/) for the web framework
- [Model Context Protocol](https://github.com/tadata-org/fastapi_mcp) for AI integration