// Hook for getting AI-powered hints
import { useState } from 'react'
import { mathService } from '../services/mathService'

export function useHint() {
  const [hint, setHint] = useState<string | null>(null)
  const [isLoading, setIsLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)

  const getHint = async (problem: string, currentStep?: string) => {
    setIsLoading(true)
    setError(null)
    
    try {
      const hintText = await mathService.getHint(problem, currentStep)
      setHint(hintText)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to get hint')
      setHint(null)
    } finally {
      setIsLoading(false)
    }
  }

  return {
    hint,
    isLoading,
    error,
    getHint
  }
}