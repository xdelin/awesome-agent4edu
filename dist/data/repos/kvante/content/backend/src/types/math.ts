// Type definitions for math-related data structures
export interface MathProblem {
  id: string;
  text: string;
  context?: string;
  timestamp: Date;
}

export interface MathSolution {
  problemId: string;
  steps: string[];
  timestamp: Date;
}

export interface MathHint {
  problem: string;
  hint: string;
  timestamp: Date;
}

export interface StepCheck {
  isCorrect: boolean;
  feedback: string;
  suggestions: string[];
}