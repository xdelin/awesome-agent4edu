// Controller for handling math problem solving requests
import { Request, Response } from 'express';
import { AIService } from '../services/AIService';
import { MathProblem, MathHint, StepCheck } from '../types/math';

export class MathController {
  private aiService: AIService;

  constructor() {
    this.aiService = new AIService();
  }

  solveProblem = async (req: Request, res: Response): Promise<void> => {
    try {
      const { problem, context }: { problem: string; context?: string } = req.body;
      
      if (!problem) {
        res.status(400).json({ error: 'Problem text is required' });
        return;
      }

      const mathProblem: MathProblem = {
        id: Date.now().toString(),
        text: problem,
        context: context || '',
        timestamp: new Date()
      };

      const solution = await this.aiService.solveMathProblem(mathProblem);
      res.json(solution);
    } catch (error) {
      console.error('Error solving math problem:', error);
      res.status(500).json({ error: 'Failed to solve math problem' });
    }
  };

  getHint = async (req: Request, res: Response): Promise<void> => {
    try {
      const { problem, currentStep }: { problem: string; currentStep?: string } = req.body;
      
      if (!problem) {
        res.status(400).json({ error: 'Problem text is required' });
        return;
      }

      const hint = await this.aiService.getHint(problem, currentStep);
      res.json({ hint });
    } catch (error) {
      console.error('Error getting hint:', error);
      res.status(500).json({ error: 'Failed to get hint' });
    }
  };

  checkStep = async (req: Request, res: Response): Promise<void> => {
    try {
      const { problem, step }: { problem: string; step: string } = req.body;
      
      if (!problem || !step) {
        res.status(400).json({ error: 'Problem and step are required' });
        return;
      }

      const stepCheck = await this.aiService.checkStep(problem, step);
      res.json(stepCheck);
    } catch (error) {
      console.error('Error checking step:', error);
      res.status(500).json({ error: 'Failed to check step' });
    }
  };
}