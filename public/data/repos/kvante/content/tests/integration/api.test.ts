// Integration tests for API endpoints
import request from 'supertest'
import express from 'express'
import { apiRoutes } from '../../backend/src/routes'

const app = express()
app.use(express.json())
app.use('/api', apiRoutes)

describe('API Integration Tests', () => {
  describe('POST /api/math/solve', () => {
    it('should solve a math problem successfully', async () => {
      const problem = {
        problem: 'Solve for x: 2x + 5 = 13',
        context: 'Basic algebra'
      }

      const response = await request(app)
        .post('/api/math/solve')
        .send(problem)
        .expect('Content-Type', /json/)

      // Note: This will fail until backend is fully implemented
      // but shows the test structure
      expect(response.status).toBe(200)
      expect(response.body).toHaveProperty('problemId')
      expect(response.body).toHaveProperty('steps')
      expect(Array.isArray(response.body.steps)).toBe(true)
    }, 10000) // Longer timeout for AI processing

    it('should return 400 for missing problem text', async () => {
      const response = await request(app)
        .post('/api/math/solve')
        .send({})

      expect(response.status).toBe(400)
      expect(response.body).toHaveProperty('error')
    })
  })

  describe('POST /api/math/hint', () => {
    it('should provide a helpful hint', async () => {
      const hintRequest = {
        problem: 'Solve for x: 2x + 5 = 13',
        currentStep: '2x + 5 = 13'
      }

      const response = await request(app)
        .post('/api/math/hint')
        .send(hintRequest)
        .expect('Content-Type', /json/)

      expect(response.status).toBe(200)
      expect(response.body).toHaveProperty('hint')
      expect(typeof response.body.hint).toBe('string')
    })

    it('should return 400 for missing problem', async () => {
      const response = await request(app)
        .post('/api/math/hint')
        .send({ currentStep: 'some step' })

      expect(response.status).toBe(400)
      expect(response.body).toHaveProperty('error')
    })
  })

  describe('POST /api/math/check-step', () => {
    it('should evaluate a student step correctly', async () => {
      const stepCheck = {
        problem: 'Solve for x: 2x + 5 = 13',
        step: '2x = 8'
      }

      const response = await request(app)
        .post('/api/math/check-step')
        .send(stepCheck)
        .expect('Content-Type', /json/)

      expect(response.status).toBe(200)
      expect(response.body).toHaveProperty('isCorrect')
      expect(response.body).toHaveProperty('feedback')
      expect(response.body).toHaveProperty('suggestions')
      expect(typeof response.body.isCorrect).toBe('boolean')
    })
  })

  describe('POST /api/feedback', () => {
    it('should submit feedback successfully', async () => {
      const feedback = {
        rating: 5,
        comment: 'Great explanation!',
        sessionId: 'test-session-123',
        problemId: 'test-problem-456'
      }

      const response = await request(app)
        .post('/api/feedback')
        .send(feedback)
        .expect('Content-Type', /json/)

      expect(response.status).toBe(200)
      expect(response.body).toHaveProperty('message')
    })

    it('should return 400 for missing rating', async () => {
      const response = await request(app)
        .post('/api/feedback')
        .send({ comment: 'No rating provided' })

      expect(response.status).toBe(400)
      expect(response.body).toHaveProperty('error')
    })
  })

  describe('GET /api/feedback', () => {
    it('should retrieve feedback with filters', async () => {
      const response = await request(app)
        .get('/api/feedback?sessionId=test-session-123')
        .expect('Content-Type', /json/)

      expect(response.status).toBe(200)
      expect(Array.isArray(response.body)).toBe(true)
    })
  })
})