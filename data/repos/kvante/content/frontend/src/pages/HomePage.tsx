// Home page component with introduction and features
import { Link } from 'react-router-dom'
import { Camera, BookOpen, MessageSquare } from 'lucide-react'

export function HomePage() {
  return (
    <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-12">
      <div className="text-center">
        <h1 className="text-4xl font-bold text-gray-900 sm:text-6xl">
          Learn Math Step by Step
        </h1>
        <p className="mt-6 text-lg text-gray-600 max-w-3xl mx-auto">
          Kvante helps students solve handwritten math problems with AI-powered step-by-step guidance. 
          Take a photo of your work and get personalized hints without giving away the answer.
        </p>
        <div className="mt-10">
          <Link
            to="/solve"
            className="bg-primary-600 text-white px-8 py-3 rounded-lg text-lg font-medium hover:bg-primary-700 transition-colors"
          >
            Start Solving Problems
          </Link>
        </div>
      </div>

      <div className="mt-20">
        <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
          <div className="text-center">
            <div className="flex justify-center mb-4">
              <Camera className="h-12 w-12 text-primary-600" />
            </div>
            <h3 className="text-xl font-semibold text-gray-900">Scan Your Work</h3>
            <p className="mt-2 text-gray-600">
              Take a photo of your handwritten math problems and let our OCR technology extract the text.
            </p>
          </div>

          <div className="text-center">
            <div className="flex justify-center mb-4">
              <BookOpen className="h-12 w-12 text-primary-600" />
            </div>
            <h3 className="text-xl font-semibold text-gray-900">Step-by-Step Guidance</h3>
            <p className="mt-2 text-gray-600">
              Get AI-powered hints and explanations that help you understand the process, not just the answer.
            </p>
          </div>

          <div className="text-center">
            <div className="flex justify-center mb-4">
              <MessageSquare className="h-12 w-12 text-primary-600" />
            </div>
            <h3 className="text-xl font-semibold text-gray-900">Learn and Improve</h3>
            <p className="mt-2 text-gray-600">
              Provide feedback to help improve the learning experience for you and other students.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}