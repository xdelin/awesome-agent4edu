// Service for math problem solving API calls
import axios from 'axios'
import { MathSolution } from '../types/math'

const api = axios.create({
  baseURL: '/api/math',
  timeout: 30000,
})

export const mathService = {
  async solveProblem(problem: string, context?: string): Promise<MathSolution> {
    const response = await api.post('/solve', { problem, context })
    return response.data
  },

  async getHint(problem: string, currentStep?: string): Promise<string> {
    const response = await api.post('/hint', { problem, currentStep })
    return response.data.hint
  },

  async checkStep(problem: string, step: string) {
    const response = await api.post('/check-step', { problem, step })
    return response.data
  }
}