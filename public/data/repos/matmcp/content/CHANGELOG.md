# Changelog

All notable changes to the Materials MCP project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial project setup with Poetry dependency management
- Basic FastAPI server implementation
- MCP protocol endpoint with `initialize` method
- GNoME OPTIMADE API exploration script
- Basic project documentation
- Git repository setup with proper `.gitignore`
- JSON-RPC 2.0 protocol support
- Pydantic models for request/response validation
- API documentation via FastAPI's Swagger UI and ReDoc
- Progress summary and changelog documentation

### Changed
- Updated Python version requirement to >=3.10,<4.0
- Disabled package mode in `pyproject.toml` for development

### Fixed
- Corrected remote repository name case sensitivity issue

### Security
- Added `.gitignore` to prevent sensitive files from being committed

## [0.1.0] - 2024-03-19

### Added
- Initial project structure
- Basic README with project description
- Core dependencies:
  - `optimade>=1.2.4`
  - `requests>=2.31.0`
  - `fastapi>=0.110.0`
  - `uvicorn>=0.27.1`
  - `pydantic>=2.6.3`

[Unreleased]: https://github.com/ZuchGuillotine/MatMCP/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/ZuchGuillotine/MatMCP/releases/tag/v0.1.0 