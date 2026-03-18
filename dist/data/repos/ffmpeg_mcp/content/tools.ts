import { z } from "zod"

export const tools = {
  speed_up: {
    name: "speed_up",
    description: "Speed up a video",
    input: {
      input_file: z.string().describe("Path to input file"),
      output_file: z
        .string()
        .describe(
          "Path to output file, output to the same directory if not specified",
        )
        .optional(),
      max_fps: z
        .number()
        .min(1)
        .max(60)
        .default(30)
        .describe("Max FPS for the output file"),
      speed_factor: z
        .number()
        .min(0.1)
        .max(10)
        .default(2)
        .describe("Speed factor for the output file, default to 2x sped up"),
    },
  },

  extract_audio: {
    name: "extract_audio",
    description: "Extract audio as mp3 from a video",
    input: {
      input_file: z.string().describe("Path to input file"),
      output_file: z
        .string()
        .describe(
          "Path to output file, output to the same directory if not specified",
        )
        .optional(),
    },
  },
}
