// Hook for OCR text extraction functionality
import { useState } from 'react'
import { ocrService } from '../services/ocrService'

export function useOCR() {
  const [isLoading, setIsLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)

  const extractText = async (file: File): Promise<string | null> => {
    setIsLoading(true)
    setError(null)
    
    try {
      const text = await ocrService.extractText(file)
      return text
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to extract text')
      return null
    } finally {
      setIsLoading(false)
    }
  }

  return {
    extractText,
    isLoading,
    error
  }
}