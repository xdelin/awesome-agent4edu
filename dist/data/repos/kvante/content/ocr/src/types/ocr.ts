// Type definitions for OCR processing
export interface OCRConfig {
  language?: string
  oem?: number // OCR Engine Mode
  psm?: number // Page Segmentation Mode
}

export interface ProcessingOptions {
  preprocessImage?: boolean
  postProcessText?: boolean
  grayscale?: boolean
  enhanceContrast?: boolean
  sharpen?: boolean
  denoise?: boolean
  threshold?: number
  resize?: {
    width: number
    height: number
  }
  removeExtraSpaces?: boolean
  filterMathContent?: boolean
}

export interface OCRResult {
  text: string
  confidence: number
  processingTime: number
  blocks?: any[]
  words?: any[]
}

export interface OCRError {
  message: string
  code: string
  details?: any
}