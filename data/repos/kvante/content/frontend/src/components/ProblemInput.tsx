// Component for manual text input of math problems
import { useState } from 'react'

interface ProblemInputProps {
  onSubmit: (problem: string) => void
}

export function ProblemInput({ onSubmit }: ProblemInputProps) {
  const [problem, setProblem] = useState('')

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault()
    if (problem.trim()) {
      onSubmit(problem.trim())
    }
  }

  return (
    <form onSubmit={handleSubmit} className="space-y-4">
      <textarea
        value={problem}
        onChange={(e) => setProblem(e.target.value)}
        placeholder="Enter your math problem here..."
        rows={4}
        className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-primary-500 focus:border-primary-500"
      />
      <button
        type="submit"
        disabled={!problem.trim()}
        className="w-full bg-primary-600 text-white px-4 py-2 rounded hover:bg-primary-700 disabled:opacity-50 disabled:cursor-not-allowed"
      >
        Solve Problem
      </button>
    </form>
  )
}