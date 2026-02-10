// Unit tests for OCR processing functionality
import { OCRProcessor } from '../../../ocr/src/core/OCRProcessor'
import { ProcessingOptions } from '../../../ocr/src/types/ocr'

// Mock Tesseract.js
jest.mock('tesseract.js', () => ({
  createWorker: jest.fn()
}))

describe('OCRProcessor', () => {
  let processor: OCRProcessor
  const mockImageBuffer = Buffer.from('fake-image-data')

  beforeEach(() => {
    processor = new OCRProcessor()
  })

  afterEach(() => {
    jest.clearAllMocks()
  })

  describe('processImage', () => {
    it('should extract text from image successfully', async () => {
      const mockWorker = {
        setParameters: jest.fn().mockResolvedValue(undefined),
        recognize: jest.fn().mockResolvedValue({
          data: {
            text: '2x + 5 = 13',
            confidence: 95,
            blocks: [],
            words: []
          }
        }),
        terminate: jest.fn().mockResolvedValue(undefined)
      }

      const { createWorker } = require('tesseract.js')
      createWorker.mockResolvedValue(mockWorker)

      const result = await processor.processImage(mockImageBuffer)

      expect(result.text).toBe('2x + 5 = 13')
      expect(result.confidence).toBe(95)
      expect(result.processingTime).toBeGreaterThan(0)
      expect(mockWorker.setParameters).toHaveBeenCalled()
      expect(mockWorker.recognize).toHaveBeenCalled()
      expect(mockWorker.terminate).toHaveBeenCalled()
    })

    it('should handle OCR errors gracefully', async () => {
      const { createWorker } = require('tesseract.js')
      createWorker.mockRejectedValue(new Error('OCR Failed'))

      await expect(processor.processImage(mockImageBuffer)).rejects.toThrow('OCR processing failed')
    })

    it('should apply preprocessing when enabled', async () => {
      const mockWorker = {
        setParameters: jest.fn().mockResolvedValue(undefined),
        recognize: jest.fn().mockResolvedValue({
          data: {
            text: 'Enhanced text',
            confidence: 98,
            blocks: [],
            words: []
          }
        }),
        terminate: jest.fn().mockResolvedValue(undefined)
      }

      const { createWorker } = require('tesseract.js')
      createWorker.mockResolvedValue(mockWorker)

      const options: ProcessingOptions = {
        preprocessImage: true,
        enhanceContrast: true,
        sharpen: true
      }

      const result = await processor.processImage(mockImageBuffer, options)

      expect(result.text).toBe('Enhanced text')
      expect(result.confidence).toBe(98)
    })
  })

  describe('extractMathText', () => {
    it('should optimize processing for math content', async () => {
      const mockWorker = {
        setParameters: jest.fn().mockResolvedValue(undefined),
        recognize: jest.fn().mockResolvedValue({
          data: {
            text: 'x² + 2x - 8 = 0',
            confidence: 92,
            blocks: [],
            words: []
          }
        }),
        terminate: jest.fn().mockResolvedValue(undefined)
      }

      const { createWorker } = require('tesseract.js')
      createWorker.mockResolvedValue(mockWorker)

      const result = await processor.extractMathText(mockImageBuffer)

      expect(result.text).toBe('x² + 2x - 8 = 0')
      expect(result.confidence).toBe(92)
    })
  })
})