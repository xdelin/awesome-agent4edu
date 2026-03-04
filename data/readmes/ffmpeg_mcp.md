# ffmpeg-mcp

A stdio MCP server to interact with ffmpeg for common media operations.

This project is sponsored by [ChatWise](https://chatwise.app) - All-in-one LLM Chatbot with MCP support.

## Usage

You need [ffmpeg](https://www.ffmpeg.org/) installed on your system first.

Then add the following command to your MCP client:

```
npx -y ffmpeg-mcp
```

By default it uses `ffmpeg` from your system path. You can also specify the path to `ffmpeg` by setting the `FFMPEG_PATH` environment variable, like:

```
FFMPEG_PATH=/path/to/ffmpeg
```

## Tools

View [tools.ts](./tools.ts)

More tools to come.

## License

MIT.
