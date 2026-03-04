# Performance

PDF Reader MCP is optimized for speed and efficiency.

## Benchmarks

Benchmarks run on Node.js 22, measuring operations per second.

| Operation | Ops/sec | Notes |
|-----------|---------|-------|
| Metadata only | ~5,000 | Fastest extraction mode |
| Single page text | ~5,300 | Minimal parsing |
| Full text (10 pages) | ~4,500 | Depends on content |
| With images | ~2,000 | Image encoding overhead |

## Optimization Tips

### 1. Request Only What You Need

```json
// Fast - metadata only
{
  "sources": [{ "path": "doc.pdf" }],
  "include_metadata": true,
  "include_page_count": true,
  "include_full_text": false,
  "include_images": false
}
```

### 2. Use Page Ranges

Instead of full text extraction, request specific pages:

```json
{
  "sources": [{
    "path": "doc.pdf",
    "pages": [1, 2]  // Only first two pages
  }],
  "include_full_text": false,
  "include_metadata": false,
  "include_page_count": false,
  "include_images": false
}
```

### 3. Batch Sources

Process multiple PDFs in one request for better throughput:

```json
{
  "sources": [
    { "path": "doc1.pdf" },
    { "path": "doc2.pdf" },
    { "path": "doc3.pdf" }
  ],
  "include_full_text": true,
  "include_metadata": false,
  "include_page_count": false,
  "include_images": false
}
```

### 4. Avoid Images Unless Needed

Image extraction involves encoding to PNG and base64, which adds overhead:

```json
// Slower
{ "include_images": true }

// Faster
{ "include_images": false }
```

## Concurrency

The server processes multiple sources concurrently with a default limit of 3 simultaneous operations to prevent memory exhaustion.

## File Size Limits

- Maximum file size: 100MB
- Files exceeding this limit will return an error

## Memory Usage

Memory usage scales with:
- Number of concurrent sources
- PDF complexity
- Image extraction enabled

For large PDFs or many concurrent requests, ensure adequate system memory.
