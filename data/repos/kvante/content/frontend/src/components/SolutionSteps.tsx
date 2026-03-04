// Component for displaying step-by-step solution guidance
import { useState } from 'react'
import { ChevronRight, Lightbulb } from 'lucide-react'
import { useHint } from '../hooks/useHint'

interface SolutionStepsProps {
  steps: string[]
  problemText: string
}

export function SolutionSteps({ steps, problemText }: SolutionStepsProps) {
  const [currentStep, setCurrentStep] = useState(0)
  const [showHint, setShowHint] = useState(false)
  const { hint, isLoading: isLoadingHint, getHint } = useHint()

  const handleNextStep = () => {
    if (currentStep < steps.length - 1) {
      setCurrentStep(currentStep + 1)
    }
  }

  const handleGetHint = async () => {
    const currentStepText = steps[currentStep]
    await getHint(problemText, currentStepText)
    setShowHint(true)
  }

  return (
    <div className="bg-white rounded-lg shadow-sm border p-6">
      <h2 className="text-xl font-semibold mb-6">Step-by-Step Guidance</h2>
      
      <div className="space-y-6">
        <div className="flex items-center space-x-4">
          <div className="flex-shrink-0 w-8 h-8 bg-primary-600 text-white rounded-full flex items-center justify-center font-semibold">
            {currentStep + 1}
          </div>
          <div className="flex-1">
            <p className="text-lg text-gray-800">{steps[currentStep]}</p>
          </div>
        </div>

        <div className="flex space-x-4">
          <button
            onClick={handleGetHint}
            disabled={isLoadingHint}
            className="inline-flex items-center px-4 py-2 border border-gray-300 rounded-md text-gray-700 hover:bg-gray-50"
          >
            <Lightbulb className="h-4 w-4 mr-2" />
            {isLoadingHint ? 'Getting hint...' : 'Get Hint'}
          </button>
          
          {currentStep < steps.length - 1 && (
            <button
              onClick={handleNextStep}
              className="inline-flex items-center px-4 py-2 bg-primary-600 text-white rounded-md hover:bg-primary-700"
            >
              Next Step
              <ChevronRight className="h-4 w-4 ml-2" />
            </button>
          )}
        </div>

        {showHint && hint && (
          <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-4">
            <h4 className="font-semibold text-yellow-800 mb-2">Hint:</h4>
            <p className="text-yellow-700">{hint}</p>
          </div>
        )}

        <div className="border-t pt-4">
          <div className="text-sm text-gray-500 mb-2">
            Step {currentStep + 1} of {steps.length}
          </div>
          <div className="w-full bg-gray-200 rounded-full h-2">
            <div 
              className="bg-primary-600 h-2 rounded-full transition-all duration-300"
              style={{ width: `${((currentStep + 1) / steps.length) * 100}%` }}
            />
          </div>
        </div>
      </div>
    </div>
  )
}