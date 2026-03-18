import { McpError, ErrorCode } from "@modelcontextprotocol/sdk/types.js";
import { AxiosInstance, AxiosError, AxiosHeaders } from "axios";
import { describe, it, expect, vi, Mock } from "vitest";

import { handleGetAuthorInfo } from "./index.js";

// Mock Axios instance
const mockAxiosInstance = {
  get: vi.fn(),
} as unknown as AxiosInstance;

describe("handleGetAuthorInfo", () => {
  it("should return author info for a valid author key", async () => {
    const mockAuthorData = {
      name: "J.R.R. Tolkien",
      key: "/authors/OL26320A",
      birth_date: "3 January 1892",
      death_date: "2 September 1973",
      bio: "British writer, poet, philologist, and academic.",
    };
    (mockAxiosInstance.get as Mock).mockResolvedValue({
      data: mockAuthorData,
    });

    const result = await handleGetAuthorInfo(
      { author_key: "OL26320A" },
      mockAxiosInstance,
    );

    expect(mockAxiosInstance.get).toHaveBeenCalledWith(
      "/authors/OL26320A.json",
    );
    expect(result.isError).toBeUndefined();
    expect(result.content).toEqual([
      {
        type: "text",
        text: JSON.stringify(mockAuthorData, null, 2),
      },
    ]);
  });

  it("should handle bio as an object", async () => {
    const mockAuthorDataWithObjectBio = {
      name: "George Orwell",
      key: "/authors/OL27346A",
      bio: {
        type: "/type/text",
        value: "English novelist, essayist, journalist and critic.",
      },
    };
    const expectedFormattedData = {
      name: "George Orwell",
      key: "/authors/OL27346A",
      bio: "English novelist, essayist, journalist and critic.",
    };
    (mockAxiosInstance.get as Mock).mockResolvedValue({
      data: mockAuthorDataWithObjectBio,
    });

    const result = await handleGetAuthorInfo(
      { author_key: "OL27346A" },
      mockAxiosInstance,
    );

    expect(result.isError).toBeUndefined();
    expect(result.content).toEqual([
      {
        type: "text",
        text: JSON.stringify(expectedFormattedData, null, 2),
      },
    ]);
  });

  it("should throw McpError for invalid author key format", async () => {
    await expect(
      handleGetAuthorInfo({ author_key: "invalid-key" }, mockAxiosInstance),
    ).rejects.toThrow(
      new McpError(
        ErrorCode.InvalidParams,
        "Invalid arguments for get_author_info: author_key: Author key must be in the format OL<number>A",
      ),
    );
  });

  it("should throw McpError for empty author key", async () => {
    await expect(
      handleGetAuthorInfo({ author_key: "" }, mockAxiosInstance),
    ).rejects.toThrow(
      new McpError(
        ErrorCode.InvalidParams,
        "Invalid arguments for get_author_info: author_key: Author key cannot be empty, author_key: Author key must be in the format OL<number>A",
      ),
    );
  });

  it("should return an error message for a 404 Not Found response", async () => {
    const authorKey = "OL00000A";
    const axiosError = new AxiosError(
      `Request failed with status code 404`,
      "404",
      undefined,
      undefined,
      {
        status: 404,
        statusText: "Not Found",
        headers: {},
        config: { headers: new AxiosHeaders() },
        data: {},
      },
    );
    (mockAxiosInstance.get as Mock).mockRejectedValue(axiosError);

    const result = await handleGetAuthorInfo(
      { author_key: authorKey },
      mockAxiosInstance,
    );

    expect(result.isError).toBe(true);
    expect(result.content).toEqual([
      {
        type: "text",
        text: `Author with key "${authorKey}" not found.`,
      },
    ]);
  });

  it("should return a generic API error message for other Axios errors", async () => {
    const authorKey = "OL12345A";
    const axiosError = new AxiosError(
      `Request failed with status code 500`,
      "500",
      undefined,
      undefined,
      {
        status: 500,
        statusText: "Internal Server Error",
        headers: {},
        config: { headers: new AxiosHeaders() },
        data: {},
      },
    );
    (mockAxiosInstance.get as Mock).mockRejectedValue(axiosError);

    const result = await handleGetAuthorInfo(
      { author_key: authorKey },
      mockAxiosInstance,
    );

    expect(result.isError).toBe(true);
    expect(result.content).toEqual([
      {
        type: "text",
        text: `Open Library API error: Internal Server Error`,
      },
    ]);
  });

  it("should return a generic error message for non-Axios errors", async () => {
    const authorKey = "OL98765A";
    const genericError = new Error("Network Error");
    (mockAxiosInstance.get as Mock).mockRejectedValue(genericError);

    const result = await handleGetAuthorInfo(
      { author_key: authorKey },
      mockAxiosInstance,
    );

    expect(result.isError).toBe(true);
    expect(result.content).toEqual([
      {
        type: "text",
        text: `Error processing request: Network Error`,
      },
    ]);
  });

  it("should return a message if API returns 200 but no data", async () => {
    const authorKey = "OL11111A";
    (mockAxiosInstance.get as Mock).mockResolvedValue({ data: null }); // Simulate no data

    const result = await handleGetAuthorInfo(
      { author_key: authorKey },
      mockAxiosInstance,
    );

    expect(result.isError).toBeUndefined();
    expect(result.content).toEqual([
      {
        type: "text",
        text: `No data found for author key: "${authorKey}"`,
      },
    ]);
  });
});
