// Text post-processing utilities for cleaning OCR output
import { ProcessingOptions } from '../types/ocr'

export class TextPostProcessor {
  process(text: string, options: ProcessingOptions = {}): string {
    let processed = text

    processed = this.cleanWhitespace(processed)
    processed = this.fixCommonOCRErrors(processed)

    if (options.filterMathContent) {
      processed = this.enhanceMathContent(processed)
    }

    if (options.removeExtraSpaces !== false) {
      processed = this.normalizeSpaces(processed)
    }

    return processed.trim()
  }

  private cleanWhitespace(text: string): string {
    return text
      .replace(/\r\n/g, '\n')
      .replace(/\r/g, '\n')
      .replace(/\n+/g, '\n')
      .replace(/[ \t]+/g, ' ')
  }

  private fixCommonOCRErrors(text: string): string {
    const corrections: Record<string, string> = {
      '0': '0',
      'O': '0',
      'l': '1',
      'I': '1',
      '|': '1',
      'S': '5',
      'B': '8',
      'Z': '2',
      'G': '6',
      'o': '0'
    }

    let corrected = text
    for (const [wrong, right] of Object.entries(corrections)) {
      const regex = new RegExp(`\\b${wrong}\\b`, 'g')
      corrected = corrected.replace(regex, right)
    }

    return corrected
  }

  private enhanceMathContent(text: string): string {
    return text
      .replace(/\s*\+\s*/g, ' + ')
      .replace(/\s*-\s*/g, ' - ')
      .replace(/\s*\*\s*/g, ' × ')
      .replace(/\s*\/\s*/g, ' ÷ ')
      .replace(/\s*=\s*/g, ' = ')
      .replace(/\s*<\s*/g, ' < ')
      .replace(/\s*>\s*/g, ' > ')
      .replace(/\(\s*/g, '(')
      .replace(/\s*\)/g, ')')
      .replace(/\[\s*/g, '[')
      .replace(/\s*\]/g, ']')
      .replace(/{\s*/g, '{')
      .replace(/\s*}/g, '}')
  }

  private normalizeSpaces(text: string): string {
    return text.replace(/\s+/g, ' ')
  }
}