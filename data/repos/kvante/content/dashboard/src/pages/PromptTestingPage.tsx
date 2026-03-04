// Prompt testing and evaluation interface
import { useState } from 'react'
import { Play, Save, RotateCcw, FileText } from 'lucide-react'

const testProblems = [
  {
    id: '1',
    category: 'Algebra',
    problem: 'Solve for x: 2x + 5 = 13',
    difficulty: 'Basic'
  },
  {
    id: '2',
    category: 'Geometry',
    problem: 'Find the area of a circle with radius 5 units',
    difficulty: 'Intermediate'
  },
  {
    id: '3',
    category: 'Calculus',
    problem: 'Find the derivative of f(x) = x² + 3x - 2',
    difficulty: 'Advanced'
  }
]

const promptTemplates = [
  {
    id: 'step-by-step-v1',
    name: 'Step-by-Step Solver v1',
    description: 'Basic step-by-step guidance',
    template: `You are a math tutor. Guide the student through solving: {problem}

Provide step-by-step guidance without giving the final answer.`
  },
  {
    id: 'step-by-step-v2',
    name: 'Step-by-Step Solver v2',
    description: 'Enhanced with questioning',
    template: `You are a helpful math tutor. Help the student solve: {problem}

Use the Socratic method - ask guiding questions and provide hints, but never give the final answer directly.`
  }
]

export function PromptTestingPage() {
  const [selectedProblem, setSelectedProblem] = useState(testProblems[0])
  const [selectedPrompt, setSelectedPrompt] = useState(promptTemplates[0])
  const [customPrompt, setCustomPrompt] = useState('')
  const [result, setResult] = useState('')
  const [isLoading, setIsLoading] = useState(false)
  const [evaluation, setEvaluation] = useState({
    accuracy: 0,
    helpfulness: 0,
    clarity: 0,
    engagement: 0
  })

  const handleRunTest = async () => {
    setIsLoading(true)
    
    // Simulate API call
    setTimeout(() => {
      setResult(`Mock response for problem: "${selectedProblem.problem}"

I see you're working with a linear equation. Let's think about this step by step.

What's the first thing you notice about this equation? We have a variable term (2x) and a constant term (5) on the left side, and just a number (13) on the right side.

What operation could we perform on both sides to start isolating the x term?

Think about the order of operations and what we need to "undo" to get x by itself.`)
      
      setIsLoading(false)
    }, 2000)
  }

  const handleSaveTest = () => {
    console.log('Saving test result...', {
      problem: selectedProblem,
      prompt: selectedPrompt,
      result,
      evaluation
    })
  }

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <h1 className="text-3xl font-bold text-gray-900">Prompt Testing</h1>
        <div className="flex space-x-3">
          <button
            onClick={handleSaveTest}
            disabled={!result}
            className="inline-flex items-center px-4 py-2 border border-gray-300 rounded-md text-sm font-medium text-gray-700 bg-white hover:bg-gray-50 disabled:opacity-50"
          >
            <Save className="h-4 w-4 mr-2" />
            Save Test
          </button>
          <button
            onClick={() => {
              setResult('')
              setEvaluation({ accuracy: 0, helpfulness: 0, clarity: 0, engagement: 0 })
            }}
            className="inline-flex items-center px-4 py-2 border border-gray-300 rounded-md text-sm font-medium text-gray-700 bg-white hover:bg-gray-50"
          >
            <RotateCcw className="h-4 w-4 mr-2" />
            Reset
          </button>
        </div>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Input Panel */}
        <div className="space-y-6">
          <div className="bg-white rounded-lg border p-6">
            <h3 className="text-lg font-semibold text-gray-900 mb-4">Test Problem</h3>
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Select Problem
                </label>
                <select
                  value={selectedProblem.id}
                  onChange={(e) => setSelectedProblem(testProblems.find(p => p.id === e.target.value)!)}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-primary-500"
                >
                  {testProblems.map((problem) => (
                    <option key={problem.id} value={problem.id}>
                      {problem.category} - {problem.problem.substring(0, 50)}...
                    </option>
                  ))}
                </select>
              </div>
              
              <div className="bg-gray-50 p-4 rounded">
                <p className="text-sm text-gray-600 mb-2">
                  <strong>Category:</strong> {selectedProblem.category} | 
                  <strong> Difficulty:</strong> {selectedProblem.difficulty}
                </p>
                <p className="text-gray-900">{selectedProblem.problem}</p>
              </div>
            </div>
          </div>

          <div className="bg-white rounded-lg border p-6">
            <h3 className="text-lg font-semibold text-gray-900 mb-4">Prompt Template</h3>
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Select Template
                </label>
                <select
                  value={selectedPrompt.id}
                  onChange={(e) => setSelectedPrompt(promptTemplates.find(p => p.id === e.target.value)!)}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-primary-500"
                >
                  {promptTemplates.map((template) => (
                    <option key={template.id} value={template.id}>
                      {template.name}
                    </option>
                  ))}
                </select>
              </div>
              
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-2">
                  Custom Prompt (optional)
                </label>
                <textarea
                  value={customPrompt}
                  onChange={(e) => setCustomPrompt(e.target.value)}
                  rows={6}
                  placeholder="Enter custom prompt or modify the template..."
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-primary-500"
                />
              </div>
              
              <button
                onClick={handleRunTest}
                disabled={isLoading}
                className="w-full bg-primary-600 text-white px-4 py-2 rounded-md hover:bg-primary-700 disabled:opacity-50 flex items-center justify-center"
              >
                {isLoading ? (
                  <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-white mr-2"></div>
                ) : (
                  <Play className="h-4 w-4 mr-2" />
                )}
                {isLoading ? 'Running Test...' : 'Run Test'}
              </button>
            </div>
          </div>
        </div>

        {/* Results Panel */}
        <div className="space-y-6">
          {result && (
            <>
              <div className="bg-white rounded-lg border p-6">
                <h3 className="text-lg font-semibold text-gray-900 mb-4 flex items-center">
                  <FileText className="h-5 w-5 mr-2" />
                  AI Response
                </h3>
                <div className="bg-gray-50 p-4 rounded max-h-64 overflow-y-auto">
                  <pre className="whitespace-pre-wrap text-sm text-gray-900">{result}</pre>
                </div>
              </div>

              <div className="bg-white rounded-lg border p-6">
                <h3 className="text-lg font-semibold text-gray-900 mb-4">Evaluation</h3>
                <div className="space-y-4">
                  {Object.entries(evaluation).map(([metric, value]) => (
                    <div key={metric}>
                      <div className="flex justify-between items-center mb-2">
                        <label className="text-sm font-medium text-gray-700 capitalize">
                          {metric}
                        </label>
                        <span className="text-sm text-gray-500">{value}/5</span>
                      </div>
                      <input
                        type="range"
                        min="0"
                        max="5"
                        step="1"
                        value={value}
                        onChange={(e) => setEvaluation(prev => ({
                          ...prev,
                          [metric]: parseInt(e.target.value)
                        }))}
                        className="w-full"
                      />
                    </div>
                  ))}
                  
                  <div className="mt-4 pt-4 border-t">
                    <div className="text-sm font-medium text-gray-700">
                      Overall Score: {(Object.values(evaluation).reduce((sum, val) => sum + val, 0) / 4).toFixed(1)}/5
                    </div>
                  </div>
                </div>
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  )
}