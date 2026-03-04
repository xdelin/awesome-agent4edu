// Hook for submitting user feedback
import { useState } from 'react'
import { feedbackService } from '../services/feedbackService'
import { FeedbackSubmission } from '../types/feedback'

export function useFeedback() {
  const [isLoading, setIsLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)

  const submitFeedback = async (feedback: FeedbackSubmission) => {
    setIsLoading(true)
    setError(null)
    
    try {
      await feedbackService.submitFeedback(feedback)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to submit feedback')
      throw err
    } finally {
      setIsLoading(false)
    }
  }

  return {
    submitFeedback,
    isLoading,
    error
  }
}