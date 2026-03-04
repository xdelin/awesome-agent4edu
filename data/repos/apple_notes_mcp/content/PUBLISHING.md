# Publishing Guide

This document explains how to publish the Apple Notes MCP server to PyPI.

## Prerequisites

1. **PyPI Account**: Create an account on [PyPI](https://pypi.org/account/register/)
2. **TestPyPI Account**: Create an account on [TestPyPI](https://test.pypi.org/account/register/)
3. **API Tokens**: Generate API tokens for both PyPI and TestPyPI

## Publishing Steps

### 1. Build the Package

```bash
# Install build tools
uv add --dev build twine

# Build the package
uv run python -m build
```

### 2. Test on TestPyPI

```bash
# Upload to TestPyPI
uv run twine upload --repository testpypi dist/*

# Test installation
uvx mcp-apple-notes@latest --index-url https://test.pypi.org/simple/
```

### 3. Publish to PyPI

```bash
# Upload to PyPI
uv run twine upload dist/*
```

### 4. Verify Installation

```bash
# Test the published package
uvx mcp-apple-notes@latest
```

## Version Management

### Update Version

1. Update version in `pyproject.toml`:
   ```toml
   version = "0.1.1"
   ```

2. Create a git tag:
   ```bash
   git tag v0.1.1
   git push origin v0.1.1
   ```

### Release Notes

Update the README.md with:
- New features
- Bug fixes
- Breaking changes
- Migration guide (if needed)

## Automated Publishing

### GitHub Actions

Create `.github/workflows/publish.yml`:

```yaml
name: Publish to PyPI

on:
  push:
    tags:
      - 'v*'

jobs:
  publish:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      - name: Install uv
        run: |
          curl -LsSf https://astral.sh/uv/install.sh | sh
          echo "$HOME/.local/bin" >> $GITHUB_PATH
      - name: Install dependencies
        run: uv sync
      - name: Build package
        run: uv run python -m build
      - name: Publish to PyPI
        run: uv run twine upload dist/*
        env:
          TWINE_USERNAME: __token__
          TWINE_PASSWORD: ${{ secrets.PYPI_API_TOKEN }}
```

## Troubleshooting

### Common Issues

1. **Authentication Failed**
   - Verify API token is correct
   - Ensure token has upload permissions

2. **Package Already Exists**
   - Increment version number
   - Remove old distribution files: `rm -rf dist/ build/`

3. **Build Errors**
   - Check `pyproject.toml` syntax
   - Verify all dependencies are listed

### Security

- Never commit API tokens to version control
- Use environment variables or secrets
- Rotate tokens regularly

## References

- [PyPI Packaging Guide](https://packaging.python.org/tutorials/packaging-projects/)
- [Twine Documentation](https://twine.readthedocs.io/)
- [GitHub Actions for Python](https://docs.github.com/en/actions/automating-builds-and-tests/building-and-testing-python)
