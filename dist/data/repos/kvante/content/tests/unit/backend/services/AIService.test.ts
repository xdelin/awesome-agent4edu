// Unit tests for AI service functionality
import { AIService } from '../../../../backend/src/services/AIService'
import { MathProblem } from '../../../../backend/src/types/math'

// Mock Anthropic SDK
jest.mock('@anthropic-ai/sdk', () => {
  return {
    __esModule: true,
    default: jest.fn().mockImplementation(() => ({
      messages: {
        create: jest.fn()
      }
    }))
  }
})

describe('AIService', () => {
  let aiService: AIService
  const mockProblem: MathProblem = {
    id: 'test-1',
    text: 'Solve for x: 2x + 5 = 13',
    context: 'Basic algebra problem',
    timestamp: new Date()
  }

  beforeEach(() => {
    aiService = new AIService()
  })

  afterEach(() => {
    jest.clearAllMocks()
  })

  describe('solveMathProblem', () => {
    it('should return step-by-step solution', async () => {
      const mockResponse = {
        content: [{
          type: 'text',
          text: '1. Subtract 5 from both sides\n2. Divide both sides by 2\n3. Check your answer'
        }]
      }

      const anthropicMock = require('@anthropic-ai/sdk').default
      const mockCreate = jest.fn().mockResolvedValue(mockResponse)
      anthropicMock.mockImplementation(() => ({
        messages: { create: mockCreate }
      }))

      const result = await aiService.solveMathProblem(mockProblem)

      expect(result).toBeDefined()
      expect(result.problemId).toBe(mockProblem.id)
      expect(result.steps).toHaveLength(3)
      expect(result.steps[0]).toContain('Subtract 5')
      expect(mockCreate).toHaveBeenCalledTimes(1)
    })

    it('should handle API errors gracefully', async () => {
      const anthropicMock = require('@anthropic-ai/sdk').default
      const mockCreate = jest.fn().mockRejectedValue(new Error('API Error'))
      anthropicMock.mockImplementation(() => ({
        messages: { create: mockCreate }
      }))

      await expect(aiService.solveMathProblem(mockProblem)).rejects.toThrow('Failed to solve math problem')
    })
  })

  describe('getHint', () => {
    it('should provide helpful hint without giving answer', async () => {
      const mockResponse = {
        content: [{
          type: 'text',
          text: 'Think about what operation would help you isolate the x term. What is being done to x in this equation?'
        }]
      }

      const anthropicMock = require('@anthropic-ai/sdk').default
      const mockCreate = jest.fn().mockResolvedValue(mockResponse)
      anthropicMock.mockImplementation(() => ({
        messages: { create: mockCreate }
      }))

      const hint = await aiService.getHint(mockProblem.text, '2x + 5 = 13')

      expect(hint).toBeDefined()
      expect(hint).toContain('isolate')
      expect(hint).not.toContain('x = 4') // Should not give away answer
      expect(mockCreate).toHaveBeenCalledTimes(1)
    })
  })

  describe('checkStep', () => {
    it('should correctly evaluate a valid step', async () => {
      const mockResponse = {
        content: [{
          type: 'text',
          text: 'Correct! You properly subtracted 5 from both sides to get 2x = 8. Good work!'
        }]
      }

      const anthropicMock = require('@anthropic-ai/sdk').default
      const mockCreate = jest.fn().mockResolvedValue(mockResponse)
      anthropicMock.mockImplementation(() => ({
        messages: { create: mockCreate }
      }))

      const result = await aiService.checkStep(mockProblem.text, '2x = 8')

      expect(result.isCorrect).toBe(true)
      expect(result.feedback).toContain('Correct')
      expect(result.suggestions).toBeDefined()
    })

    it('should provide constructive feedback for incorrect step', async () => {
      const mockResponse = {
        content: [{
          type: 'text',
          text: 'This step needs attention. When you divide by 2, make sure to divide each term separately. Try again!'
        }]
      }

      const anthropicMock = require('@anthropic-ai/sdk').default
      const mockCreate = jest.fn().mockResolvedValue(mockResponse)
      anthropicMock.mockImplementation(() => ({
        messages: { create: mockCreate }
      }))

      const result = await aiService.checkStep(mockProblem.text, 'x + 2.5 = 6.5')

      expect(result.isCorrect).toBe(false)
      expect(result.feedback).toContain('needs attention')
      expect(result.suggestions).toBeDefined()
    })
  })
})