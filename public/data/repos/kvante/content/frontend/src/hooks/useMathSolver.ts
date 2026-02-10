// Hook for math problem solving functionality
import { useState } from 'react'
import { mathService } from '../services/mathService'
import { MathSolution } from '../types/math'

export function useMathSolver() {
  const [solution, setSolution] = useState<MathSolution | null>(null)
  const [isLoading, setIsLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)

  const solveProblem = async (problem: string, context?: string) => {
    setIsLoading(true)
    setError(null)
    
    try {
      const result = await mathService.solveProblem(problem, context)
      setSolution(result)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to solve problem')
      setSolution(null)
    } finally {
      setIsLoading(false)
    }
  }

  const reset = () => {
    setSolution(null)
    setError(null)
  }

  return {
    solution,
    isLoading,
    error,
    solveProblem,
    reset
  }
}