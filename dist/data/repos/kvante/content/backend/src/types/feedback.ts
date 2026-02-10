// Type definitions for user feedback data structures
export interface Feedback {
  id: string;
  rating: number;
  comment?: string;
  sessionId?: string;
  problemId?: string;
  timestamp: Date;
}