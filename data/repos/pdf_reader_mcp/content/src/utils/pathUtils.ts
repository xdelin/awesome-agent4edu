import path from 'node:path';
import { ErrorCode, PdfError } from './errors.js';

// Use the server's current working directory as the project root.
// This relies on the process launching the server to set the CWD correctly.
export const PROJECT_ROOT = process.cwd();

/**
 * Resolves a user-provided path, accepting both absolute and relative paths.
 * Relative paths are resolved against the current working directory (PROJECT_ROOT).
 * @param userPath The path provided by the user (absolute or relative).
 * @returns The resolved absolute path.
 * @throws {PdfError} If path is invalid.
 */
export const resolvePath = (userPath: string): string => {
  if (typeof userPath !== 'string') {
    throw new PdfError(ErrorCode.InvalidParams, 'Path must be a string.');
  }

  const normalizedUserPath = path.normalize(userPath);

  // Resolve the path (absolute paths stay as-is, relative paths resolve against PROJECT_ROOT)
  return path.isAbsolute(normalizedUserPath)
    ? normalizedUserPath
    : path.resolve(PROJECT_ROOT, normalizedUserPath);
};
