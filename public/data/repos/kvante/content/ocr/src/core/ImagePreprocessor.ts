// Image preprocessing utilities for better OCR accuracy
import sharp from 'sharp'
import { ProcessingOptions } from '../types/ocr'

export class ImagePreprocessor {
  async preprocess(imageBuffer: Buffer, options: ProcessingOptions = {}): Promise<Buffer> {
    try {
      let image = sharp(imageBuffer)

      if (options.resize) {
        image = image.resize(options.resize.width, options.resize.height, {
          fit: 'inside',
          withoutEnlargement: true
        })
      }

      if (options.enhanceContrast) {
        image = image.normalise()
      }

      if (options.sharpen) {
        image = image.sharpen()
      }

      if (options.denoise) {
        image = image.median(3)
      }

      if (options.grayscale !== false) {
        image = image.grayscale()
      }

      if (options.threshold) {
        image = image.threshold(options.threshold)
      }

      return await image.png().toBuffer()
    } catch (error) {
      console.error('Image preprocessing failed:', error)
      return imageBuffer
    }
  }

  async optimizeForMath(imageBuffer: Buffer): Promise<Buffer> {
    return this.preprocess(imageBuffer, {
      grayscale: true,
      enhanceContrast: true,
      sharpen: true,
      denoise: true,
      resize: { width: 1200, height: 1200 }
    })
  }
}