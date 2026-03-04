import { McpServer } from "@modelcontextprotocol/sdk/server/mcp.js"
import { StdioServerTransport } from "@modelcontextprotocol/sdk/server/stdio.js"
import { x } from "tinyexec"
import path from "node:path"
import { tools } from "./tools"

const ffmpeg = process.env.FFMPEG_PATH || "ffmpeg"

const getOutputFilePath = (input: string, ext: string) => {
  const input_dir = path.dirname(input)
  const input_filename = path.basename(input, path.extname(input))
  return path.join(input_dir, `${input_filename}_output${ext}`)
}

const server = new McpServer(
  {
    name: "ffmpeg-mcp",
    version: "1.0.0",
  },
  {
    capabilities: {
      logging: {},
    },
  },
)

server.tool(
  tools.speed_up.name,
  tools.speed_up.description,
  tools.speed_up.input,
  async (args) => {
    const output_file =
      args.output_file ||
      getOutputFilePath(args.input_file, path.extname(args.input_file))

    const result = x(ffmpeg, [
      `-i`,
      args.input_file,
      `-filter:v`,
      `setpts=${1 / args.speed_factor}*PTS,fps=fps=${args.max_fps}`,
      `-af`,
      `atempo=${args.speed_factor}`,
      `-y`,
      output_file,
    ])

    let output = ""

    for await (const line of result) {
      output += line
      await server.server.sendLoggingMessage({
        level: "debug",
        data: {
          line,
        },
      })
    }

    if (result.exitCode !== 0) {
      return {
        content: [
          {
            type: "text",
            text: `Error: ${output}`,
          },
        ],
        isError: true,
      }
    }

    return {
      content: [
        {
          type: "text",
          text: `Successfully sped up video to: ${output_file}`,
        },
      ],
    }
  },
)

server.tool(
  tools.extract_audio.name,
  tools.extract_audio.description,
  tools.extract_audio.input,
  async (args) => {
    const output_file =
      args.output_file || getOutputFilePath(args.input_file, ".mp3")

    const result = await x(ffmpeg, [
      `-i`,
      args.input_file,
      `-vn`,
      `-acodec`,
      `mp3`,
      `-y`,
      output_file,
    ])

    if (result.exitCode !== 0) {
      return {
        content: [
          {
            type: "text",
            text: `Error: ${result.stderr}`,
          },
        ],
        isError: true,
      }
    }

    return {
      content: [
        {
          type: "text",
          text: `Successfully extracted audio to: ${output_file}`,
        },
      ],
    }
  },
)

const transport = new StdioServerTransport()
await server.connect(transport)
