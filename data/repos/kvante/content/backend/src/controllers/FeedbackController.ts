// Controller for handling user feedback collection
import { Request, Response } from 'express';
import { FeedbackService } from '../services/FeedbackService';
import { Feedback } from '../types/feedback';

export class FeedbackController {
  private feedbackService: FeedbackService;

  constructor() {
    this.feedbackService = new FeedbackService();
  }

  submitFeedback = async (req: Request, res: Response): Promise<void> => {
    try {
      const { rating, comment, sessionId, problemId }: Feedback = req.body;
      
      if (!rating) {
        res.status(400).json({ error: 'Rating is required' });
        return;
      }

      const feedback: Feedback = {
        id: Date.now().toString(),
        rating,
        comment: comment || '',
        sessionId: sessionId || '',
        problemId: problemId || '',
        timestamp: new Date()
      };

      await this.feedbackService.saveFeedback(feedback);
      res.json({ message: 'Feedback submitted successfully' });
    } catch (error) {
      console.error('Error submitting feedback:', error);
      res.status(500).json({ error: 'Failed to submit feedback' });
    }
  };

  getFeedback = async (req: Request, res: Response): Promise<void> => {
    try {
      const { sessionId, problemId } = req.query;
      const feedback = await this.feedbackService.getFeedback(sessionId as string, problemId as string);
      res.json(feedback);
    } catch (error) {
      console.error('Error getting feedback:', error);
      res.status(500).json({ error: 'Failed to get feedback' });
    }
  };
}