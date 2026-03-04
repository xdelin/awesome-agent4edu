import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { ErrorCode, PdfError } from '../src/utils/errors.js';
import { PROJECT_ROOT, resolvePath } from '../src/utils/pathUtils.js';

describe('resolvePath Utility', () => {
  it('should resolve a valid relative path correctly', () => {
    const userPath = 'some/file.txt';
    const expectedPath = path.resolve(PROJECT_ROOT, userPath);
    expect(resolvePath(userPath)).toBe(expectedPath);
  });

  it('should resolve paths with "." correctly', () => {
    const userPath = './some/./other/file.txt';
    const expectedPath = path.resolve(PROJECT_ROOT, 'some/other/file.txt');
    expect(resolvePath(userPath)).toBe(expectedPath);
  });

  it('should resolve paths with ".." correctly', () => {
    const userPath = 'some/folder/../other/file.txt';
    const expectedPath = path.resolve(PROJECT_ROOT, 'some/other/file.txt');
    expect(resolvePath(userPath)).toBe(expectedPath);
  });

  it('should resolve relative paths that go outside PROJECT_ROOT', () => {
    const userPath = '../outside/file.txt';
    const result = resolvePath(userPath);
    expect(result).toBe(path.resolve(PROJECT_ROOT, userPath));
  });

  it('should accept absolute paths and return them normalized', () => {
    const userPath = path.resolve(PROJECT_ROOT, 'absolute/file.txt');
    expect(resolvePath(userPath)).toBe(path.normalize(userPath));
  });

  it('should accept any absolute path', () => {
    const absolutePath = path.sep === '/' ? '/etc/passwd' : 'C:\\Windows\\System32\\config.txt';
    expect(resolvePath(absolutePath)).toBe(path.normalize(absolutePath));
  });

  it('should throw PdfError for non-string input', () => {
    const userPath = 123 as unknown as string;
    expect(() => resolvePath(userPath)).toThrow(PdfError);
    expect(() => resolvePath(userPath)).toThrow('Path must be a string.');
    try {
      resolvePath(userPath);
    } catch (e) {
      expect(e).toBeInstanceOf(PdfError);
      expect((e as PdfError).code).toBe(ErrorCode.InvalidParams);
    }
  });

  it('should handle empty string input', () => {
    const userPath = '';
    const expectedPath = path.resolve(PROJECT_ROOT, '');
    expect(resolvePath(userPath)).toBe(expectedPath);
  });
});
