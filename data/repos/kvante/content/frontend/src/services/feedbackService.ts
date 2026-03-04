// Service for feedback submission API calls
import axios from 'axios'
import { FeedbackSubmission } from '../types/feedback'

const api = axios.create({
  baseURL: '/api/feedback',
  timeout: 10000,
})

export const feedbackService = {
  async submitFeedback(feedback: FeedbackSubmission): Promise<void> {
    await api.post('/', feedback)
  },

  async getFeedback(sessionId?: string, problemId?: string) {
    const params = new URLSearchParams()
    if (sessionId) params.append('sessionId', sessionId)
    if (problemId) params.append('problemId', problemId)
    
    const response = await api.get(`/?${params.toString()}`)
    return response.data
  }
}