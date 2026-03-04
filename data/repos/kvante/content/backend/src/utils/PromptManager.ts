// Utility class for managing AI prompts
import { promises as fs } from 'fs';
import path from 'path';

export class PromptManager {
  private promptsDir = path.join(process.cwd(), '..', 'prompts');

  getMathSolvingPrompt(problem: string, context?: string): string {
    const basePrompt = `You are a math tutor helping a student learn step-by-step problem solving. 
Never provide the final answer directly. Instead, guide the student through the solution process with clear, educational steps.

Problem: ${problem}
${context ? `Context: ${context}` : ''}

Provide step-by-step guidance that helps the student understand the process, but does not give away the final answer.`;

    return basePrompt;
  }

  getHintPrompt(problem: string, currentStep?: string): string {
    return `You are a math tutor. The student is working on this problem: "${problem}"
${currentStep ? `They are currently at this step: "${currentStep}"` : ''}

Provide a helpful hint to guide them to the next step, but do not solve the problem for them.`;
  }

  getStepCheckPrompt(problem: string, step: string): string {
    return `You are a math tutor checking a student's work.

Problem: ${problem}
Student's step: ${step}

Evaluate if this step is correct and provide constructive feedback. If incorrect, suggest improvements without giving away the answer.`;
  }

  async loadPromptFromFile(filename: string): Promise<string> {
    try {
      const promptPath = path.join(this.promptsDir, filename);
      return await fs.readFile(promptPath, 'utf-8');
    } catch (error) {
      console.error(`Error loading prompt from ${filename}:`, error);
      return '';
    }
  }
}