# Getting Started

Once installed, the PDF Reader MCP server provides a single tool: `read_pdf`.

## Basic Usage

### Get Metadata and Page Count

```json
{
  "sources": [{ "path": "/path/to/document.pdf" }],
  "include_full_text": false,
  "include_metadata": true,
  "include_page_count": true,
  "include_images": false
}
```

### Get Full Text

```json
{
  "sources": [{ "path": "/path/to/document.pdf" }],
  "include_full_text": true,
  "include_metadata": false,
  "include_page_count": false,
  "include_images": false
}
```

### Get Specific Pages

```json
{
  "sources": [{
    "path": "/path/to/document.pdf",
    "pages": [1, 3, 5]
  }],
  "include_full_text": false,
  "include_metadata": false,
  "include_page_count": false,
  "include_images": false
}
```

Or use page ranges:

```json
{
  "sources": [{
    "path": "/path/to/document.pdf",
    "pages": "1-5, 10, 15-20"
  }],
  "include_full_text": false,
  "include_metadata": false,
  "include_page_count": false,
  "include_images": false
}
```

### Extract Images

```json
{
  "sources": [{ "path": "/path/to/document.pdf" }],
  "include_full_text": false,
  "include_metadata": false,
  "include_page_count": false,
  "include_images": true
}
```

## Multiple Sources

Process multiple PDFs in a single request:

```json
{
  "sources": [
    { "path": "/path/to/report.pdf" },
    { "path": "/path/to/invoice.pdf" },
    { "url": "https://example.com/whitepaper.pdf" }
  ],
  "include_full_text": true,
  "include_metadata": true,
  "include_page_count": true,
  "include_images": false
}
```

## Response Format

```json
{
  "results": [
    {
      "source": "/path/to/document.pdf",
      "success": true,
      "data": {
        "num_pages": 10,
        "info": {
          "Title": "Document Title",
          "Author": "Author Name",
          "CreationDate": "D:20231201120000"
        },
        "metadata": { ... },
        "page_texts": [
          { "page": 1, "text": "Page 1 content..." },
          { "page": 2, "text": "Page 2 content..." }
        ],
        "images": [
          {
            "page": 1,
            "index": 0,
            "width": 800,
            "height": 600,
            "data": "data:image/png;base64,..."
          }
        ]
      }
    }
  ]
}
```

## Error Handling

If a source fails, it will be included in results with `success: false`:

```json
{
  "results": [
    {
      "source": "/path/to/missing.pdf",
      "success": false,
      "error": {
        "code": "FileNotFound",
        "message": "File not found: /path/to/missing.pdf"
      }
    }
  ]
}
```

Other sources in the same request will still be processed.
