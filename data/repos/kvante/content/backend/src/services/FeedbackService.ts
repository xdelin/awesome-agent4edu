// Service for managing user feedback storage and retrieval
import { promises as fs } from 'fs';
import path from 'path';
import { Feedback } from '../types/feedback';

export class FeedbackService {
  private dataDir = path.join(process.cwd(), '..', 'data');
  private feedbackFile = path.join(this.dataDir, 'feedback.json');

  async saveFeedback(feedback: Feedback): Promise<void> {
    try {
      await this.ensureDataDir();
      
      let existingFeedback: Feedback[] = [];
      try {
        const data = await fs.readFile(this.feedbackFile, 'utf-8');
        existingFeedback = JSON.parse(data);
      } catch (error) {
        // File doesn't exist yet, start with empty array
      }

      existingFeedback.push(feedback);
      await fs.writeFile(this.feedbackFile, JSON.stringify(existingFeedback, null, 2));
    } catch (error) {
      console.error('Error saving feedback:', error);
      throw new Error('Failed to save feedback');
    }
  }

  async getFeedback(sessionId?: string, problemId?: string): Promise<Feedback[]> {
    try {
      const data = await fs.readFile(this.feedbackFile, 'utf-8');
      let feedback: Feedback[] = JSON.parse(data);

      if (sessionId) {
        feedback = feedback.filter(f => f.sessionId === sessionId);
      }

      if (problemId) {
        feedback = feedback.filter(f => f.problemId === problemId);
      }

      return feedback;
    } catch (error) {
      console.error('Error getting feedback:', error);
      return [];
    }
  }

  private async ensureDataDir(): Promise<void> {
    try {
      await fs.access(this.dataDir);
    } catch {
      await fs.mkdir(this.dataDir, { recursive: true });
    }
  }
}