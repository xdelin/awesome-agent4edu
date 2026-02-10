// End-to-end tests for problem solving workflow
import { test, expect } from '@playwright/test'

test.describe('Problem Solving Flow', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto('/')
  })

  test('should navigate to problem solver from home page', async ({ page }) => {
    await expect(page.locator('h1')).toContainText('Learn Math Step by Step')
    
    await page.click('text=Start Solving Problems')
    await expect(page).toHaveURL('/solve')
    await expect(page.locator('h1')).toContainText('Math Problem Solver')
  })

  test('should allow text input for math problems', async ({ page }) => {
    await page.goto('/solve')
    
    const problemText = 'Solve for x: 2x + 5 = 13'
    await page.fill('textarea[placeholder*="Enter your math problem"]', problemText)
    
    await page.click('button:has-text("Solve Problem")')
    
    // Wait for problem to be displayed
    await expect(page.locator('text=Problem')).toBeVisible()
    await expect(page.locator(`text=${problemText}`)).toBeVisible()
  })

  test('should display step-by-step guidance after solving', async ({ page }) => {
    await page.goto('/solve')
    
    // Mock the API response for testing
    await page.route('/api/math/solve', async route => {
      await route.fulfill({
        status: 200,
        contentType: 'application/json',
        body: JSON.stringify({
          problemId: 'test-123',
          steps: [
            'Identify what we need to solve for: the value of x',
            'Look at the operations being performed on x',
            'Think about how to isolate x by undoing these operations'
          ],
          timestamp: new Date().toISOString()
        })
      })
    })
    
    await page.fill('textarea[placeholder*="Enter your math problem"]', 'Solve for x: 2x + 5 = 13')
    await page.click('button:has-text("Solve Problem")')
    
    // Check that step-by-step guidance appears
    await expect(page.locator('text=Step-by-Step Guidance')).toBeVisible()
    await expect(page.locator('text=Step 1 of 3')).toBeVisible()
    await expect(page.locator('text=Identify what we need to solve')).toBeVisible()
  })

  test('should allow progression through solution steps', async ({ page }) => {
    await page.goto('/solve')
    
    // Mock API response
    await page.route('/api/math/solve', async route => {
      await route.fulfill({
        status: 200,
        contentType: 'application/json',
        body: JSON.stringify({
          problemId: 'test-123',
          steps: [
            'Step 1: Subtract 5 from both sides',
            'Step 2: Divide both sides by 2',
            'Step 3: Check your answer'
          ],
          timestamp: new Date().toISOString()
        })
      })
    })
    
    await page.fill('textarea[placeholder*="Enter your math problem"]', 'Solve for x: 2x + 5 = 13')
    await page.click('button:has-text("Solve Problem")')
    
    // Should start at step 1
    await expect(page.locator('text=Step 1 of 3')).toBeVisible()
    await expect(page.locator('text=Step 1: Subtract 5')).toBeVisible()
    
    // Click next step
    await page.click('button:has-text("Next Step")')
    await expect(page.locator('text=Step 2 of 3')).toBeVisible()
    await expect(page.locator('text=Step 2: Divide both sides')).toBeVisible()
  })

  test('should provide hints when requested', async ({ page }) => {
    await page.goto('/solve')
    
    // Mock solve API
    await page.route('/api/math/solve', async route => {
      await route.fulfill({
        status: 200,
        contentType: 'application/json',
        body: JSON.stringify({
          problemId: 'test-123',
          steps: ['Think about what operation would help isolate x'],
          timestamp: new Date().toISOString()
        })
      })
    })
    
    // Mock hint API
    await page.route('/api/math/hint', async route => {
      await route.fulfill({
        status: 200,
        contentType: 'application/json',
        body: JSON.stringify({
          hint: 'Look at what operations are being performed on x. What is the opposite of addition?'
        })
      })
    })
    
    await page.fill('textarea[placeholder*="Enter your math problem"]', 'Solve for x: 2x + 5 = 13')
    await page.click('button:has-text("Solve Problem")')
    
    await page.click('button:has-text("Get Hint")')
    
    // Check that hint appears
    await expect(page.locator('text=Hint:')).toBeVisible()
    await expect(page.locator('text=Look at what operations')).toBeVisible()
  })

  test('should allow feedback submission', async ({ page }) => {
    await page.goto('/solve')
    
    // Mock solve API
    await page.route('/api/math/solve', async route => {
      await route.fulfill({
        status: 200,
        contentType: 'application/json',
        body: JSON.stringify({
          problemId: 'test-123',
          steps: ['Sample step guidance'],
          timestamp: new Date().toISOString()
        })
      })
    })
    
    // Mock feedback API
    await page.route('/api/feedback', async route => {
      await route.fulfill({
        status: 200,
        contentType: 'application/json',
        body: JSON.stringify({
          message: 'Feedback submitted successfully'
        })
      })
    })
    
    await page.fill('textarea[placeholder*="Enter your math problem"]', 'Test problem')
    await page.click('button:has-text("Solve Problem")')
    
    // Scroll to feedback section
    await page.locator('text=How was this solution?').scrollIntoViewIfNeeded()
    
    // Rate 5 stars
    await page.click('button[aria-label="5 stars"] svg, button:has(svg) >> nth=4')
    
    // Add comment
    await page.fill('textarea[placeholder*="Tell us what worked well"]', 'Great explanation!')
    
    // Submit feedback
    await page.click('button:has-text("Submit Feedback")')
    
    // Check success message
    await expect(page.locator('text=Thank you for your feedback!')).toBeVisible()
  })
})