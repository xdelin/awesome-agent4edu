// Main OCR processing class using Tesseract.js
import { createWorker } from 'tesseract.js'
import { ImagePreprocessor } from './ImagePreprocessor'
import { TextPostProcessor } from './TextPostProcessor'
import { OCRConfig, OCRResult, ProcessingOptions } from '../types/ocr'

export class OCRProcessor {
  private config: OCRConfig
  private imagePreprocessor: ImagePreprocessor
  private textPostProcessor: TextPostProcessor

  constructor(config: OCRConfig = {}) {
    this.config = {
      language: 'eng',
      oem: 1,
      psm: 6,
      ...config
    }
    this.imagePreprocessor = new ImagePreprocessor()
    this.textPostProcessor = new TextPostProcessor()
  }

  async processImage(
    imageBuffer: Buffer, 
    options: ProcessingOptions = {}
  ): Promise<OCRResult> {
    try {
      const startTime = Date.now()
      
      let processedImage = imageBuffer
      if (options.preprocessImage !== false) {
        processedImage = await this.imagePreprocessor.preprocess(imageBuffer, options)
      }

      const worker = await createWorker(this.config.language)
      
      await worker.setParameters({
        tessedit_ocr_engine_mode: this.config.oem,
        tessedit_pageseg_mode: this.config.psm,
      })

      const { data } = await worker.recognize(processedImage)
      await worker.terminate()

      let text = data.text
      if (options.postProcessText !== false) {
        text = this.textPostProcessor.process(text, options)
      }

      const processingTime = Date.now() - startTime

      return {
        text,
        confidence: data.confidence,
        processingTime,
        blocks: data.blocks,
        words: data.words
      }
    } catch (error) {
      console.error('OCR processing failed:', error)
      throw new Error(`OCR processing failed: ${error instanceof Error ? error.message : 'Unknown error'}`)
    }
  }

  async extractMathText(imageBuffer: Buffer): Promise<OCRResult> {
    return this.processImage(imageBuffer, {
      enhanceContrast: true,
      sharpen: true,
      postProcessText: true,
      filterMathContent: true
    })
  }
}