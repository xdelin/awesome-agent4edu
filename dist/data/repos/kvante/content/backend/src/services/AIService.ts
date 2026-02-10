// Service for AI-powered math problem solving using Claude/OpenAI
import Anthropic from '@anthropic-ai/sdk';
import OpenAI from 'openai';
import { MathProblem, MathSolution, StepCheck } from '../types/math';
import { PromptManager } from '../utils/PromptManager';

export class AIService {
  private anthropic: Anthropic;
  private openai: OpenAI;
  private promptManager: PromptManager;

  constructor() {
    this.anthropic = new Anthropic({
      apiKey: process.env.ANTHROPIC_API_KEY,
    });
    
    this.openai = new OpenAI({
      apiKey: process.env.OPENAI_API_KEY,
    });
    
    this.promptManager = new PromptManager();
  }

  async solveMathProblem(problem: MathProblem): Promise<MathSolution> {
    try {
      const prompt = this.promptManager.getMathSolvingPrompt(problem.text, problem.context);
      
      const response = await this.anthropic.messages.create({
        model: 'claude-3-sonnet-20240229',
        max_tokens: 1024,
        messages: [
          {
            role: 'user',
            content: prompt
          }
        ]
      });

      const content = response.content[0];
      const steps = this.parseSteps(content.type === 'text' ? content.text : '');

      return {
        problemId: problem.id,
        steps,
        timestamp: new Date()
      };
    } catch (error) {
      console.error('Error solving math problem with Claude:', error);
      throw new Error('Failed to solve math problem');
    }
  }

  async getHint(problem: string, currentStep?: string): Promise<string> {
    try {
      const prompt = this.promptManager.getHintPrompt(problem, currentStep);
      
      const response = await this.anthropic.messages.create({
        model: 'claude-3-sonnet-20240229',
        max_tokens: 512,
        messages: [
          {
            role: 'user',
            content: prompt
          }
        ]
      });

      const content = response.content[0];
      return content.type === 'text' ? content.text : '';
    } catch (error) {
      console.error('Error getting hint:', error);
      throw new Error('Failed to get hint');
    }
  }

  async checkStep(problem: string, step: string): Promise<StepCheck> {
    try {
      const prompt = this.promptManager.getStepCheckPrompt(problem, step);
      
      const response = await this.anthropic.messages.create({
        model: 'claude-3-sonnet-20240229',
        max_tokens: 512,
        messages: [
          {
            role: 'user',
            content: prompt
          }
        ]
      });

      const content = response.content[0];
      const result = content.type === 'text' ? content.text : '';
      
      return {
        isCorrect: result.toLowerCase().includes('correct'),
        feedback: result,
        suggestions: this.extractSuggestions(result)
      };
    } catch (error) {
      console.error('Error checking step:', error);
      throw new Error('Failed to check step');
    }
  }

  private parseSteps(response: string): string[] {
    return response.split('\n')
      .filter(line => line.trim())
      .map(line => line.replace(/^\d+\.\s*/, '').trim())
      .filter(step => step.length > 0);
  }

  private extractSuggestions(feedback: string): string[] {
    const lines = feedback.split('\n');
    const suggestions = lines.filter(line => 
      line.toLowerCase().includes('try') || 
      line.toLowerCase().includes('consider') ||
      line.toLowerCase().includes('suggestion')
    );
    return suggestions.map(s => s.trim());
  }
}