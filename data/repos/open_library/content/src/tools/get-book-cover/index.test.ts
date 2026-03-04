import { McpError, ErrorCode } from "@modelcontextprotocol/sdk/types.js";
import { describe, expect, it } from "vitest";

import { handleGetBookCover } from "./index.js";

describe("handleGetBookCover", () => {
  it("should generate correct URL with ISBN", async () => {
    const result = await handleGetBookCover({
      key: "ISBN",
      value: "0451526538",
    });
    expect(result.content[0].text).toBe(
      "https://covers.openlibrary.org/b/isbn/0451526538-L.jpg",
    );
  });

  it("should generate correct URL with OLID", async () => {
    const result = await handleGetBookCover({
      key: "OLID",
      value: "OL7353617M",
    });
    expect(result.content[0].text).toBe(
      "https://covers.openlibrary.org/b/olid/OL7353617M-L.jpg",
    );
  });

  it("should handle different size parameters", async () => {
    const smallSize = await handleGetBookCover({
      key: "ISBN",
      value: "0451526538",
      size: "S",
    });
    expect(smallSize.content[0].text).toBe(
      "https://covers.openlibrary.org/b/isbn/0451526538-S.jpg",
    );

    const mediumSize = await handleGetBookCover({
      key: "ISBN",
      value: "0451526538",
      size: "M",
    });
    expect(mediumSize.content[0].text).toBe(
      "https://covers.openlibrary.org/b/isbn/0451526538-M.jpg",
    );

    const largeSize = await handleGetBookCover({
      key: "ISBN",
      value: "0451526538",
      size: "L",
    });
    expect(largeSize.content[0].text).toBe(
      "https://covers.openlibrary.org/b/isbn/0451526538-L.jpg",
    );
  });

  it("should default to large size when size is null or omitted", async () => {
    const nullSize = await handleGetBookCover({
      key: "ISBN",
      value: "0451526538",
      size: null,
    });
    expect(nullSize.content[0].text).toBe(
      "https://covers.openlibrary.org/b/isbn/0451526538-L.jpg",
    );

    const omittedSize = await handleGetBookCover({
      key: "ISBN",
      value: "0451526538",
    });
    expect(omittedSize.content[0].text).toBe(
      "https://covers.openlibrary.org/b/isbn/0451526538-L.jpg",
    );
  });

  it("should throw error for invalid key", async () => {
    await expect(
      handleGetBookCover({
        key: "INVALID_KEY",
        value: "0451526538",
      }),
    ).rejects.toThrow(McpError);

    try {
      await handleGetBookCover({ key: "INVALID_KEY", value: "0451526538" });
    } catch (error) {
      expect(error).toBeInstanceOf(McpError);
      expect((error as McpError).code).toBe(ErrorCode.InvalidParams);
      expect((error as McpError).message).toContain(
        "Invalid arguments for get_book_cover",
      );
    }
  });

  it("should throw error for empty value", async () => {
    await expect(
      handleGetBookCover({
        key: "ISBN",
        value: "",
      }),
    ).rejects.toThrow(McpError);

    try {
      await handleGetBookCover({ key: "ISBN", value: "" });
    } catch (error) {
      expect(error).toBeInstanceOf(McpError);
      expect((error as McpError).code).toBe(ErrorCode.InvalidParams);
      expect((error as McpError).message).toContain("Value cannot be empty");
    }
  });

  it("should throw error for invalid size", async () => {
    await expect(
      handleGetBookCover({
        key: "ISBN",
        value: "0451526538",
        size: "XL",
      }),
    ).rejects.toThrow(McpError);
  });

  it("should throw error for missing required parameters", async () => {
    await expect(handleGetBookCover({})).rejects.toThrow(McpError);
    await expect(handleGetBookCover({ key: "ISBN" })).rejects.toThrow(McpError);
    await expect(handleGetBookCover({ value: "0451526538" })).rejects.toThrow(
      McpError,
    );
  });
});
