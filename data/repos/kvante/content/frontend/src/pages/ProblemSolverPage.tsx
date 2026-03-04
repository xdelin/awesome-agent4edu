// Problem solver page with OCR upload and step-by-step guidance
import { useState } from 'react'
import { ImageUpload } from '../components/ImageUpload'
import { ProblemInput } from '../components/ProblemInput'
import { SolutionSteps } from '../components/SolutionSteps'
import { FeedbackForm } from '../components/FeedbackForm'
import { useMathSolver } from '../hooks/useMathSolver'

export function ProblemSolverPage() {
  const [problem, setProblem] = useState('')
  const [sessionId] = useState(() => Date.now().toString())
  const { solution, isLoading, solveProblem } = useMathSolver()

  const handleProblemSubmit = async (problemText: string) => {
    setProblem(problemText)
    await solveProblem(problemText)
  }

  const handleImageUpload = async (extractedText: string) => {
    setProblem(extractedText)
    await solveProblem(extractedText)
  }

  return (
    <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
      <h1 className="text-3xl font-bold text-gray-900 mb-8">Math Problem Solver</h1>
      
      <div className="space-y-8">
        <div className="bg-white rounded-lg shadow-sm border p-6">
          <h2 className="text-xl font-semibold mb-4">Upload Image or Enter Problem</h2>
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            <div>
              <h3 className="text-lg font-medium mb-3">Take a Photo</h3>
              <ImageUpload onTextExtracted={handleImageUpload} />
            </div>
            <div>
              <h3 className="text-lg font-medium mb-3">Or Type Your Problem</h3>
              <ProblemInput onSubmit={handleProblemSubmit} />
            </div>
          </div>
        </div>

        {problem && (
          <div className="bg-white rounded-lg shadow-sm border p-6">
            <h2 className="text-xl font-semibold mb-4">Problem</h2>
            <p className="text-gray-700 bg-gray-50 p-4 rounded">{problem}</p>
          </div>
        )}

        {isLoading && (
          <div className="bg-white rounded-lg shadow-sm border p-6">
            <div className="flex items-center justify-center">
              <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-primary-600"></div>
              <span className="ml-3 text-gray-600">Generating step-by-step guidance...</span>
            </div>
          </div>
        )}

        {solution && (
          <>
            <SolutionSteps 
              steps={solution.steps}
              problemText={problem}
            />
            <FeedbackForm 
              sessionId={sessionId}
              problemId={solution.problemId}
            />
          </>
        )}
      </div>
    </div>
  )
}