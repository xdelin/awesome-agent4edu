// User feedback analysis and visualization page
import { PieChart, Pie, Cell, BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts'
import { Star, TrendingUp, MessageSquare, Filter } from 'lucide-react'

const mockFeedbackData = {
  ratingDistribution: [
    { rating: '5 Stars', count: 45, percentage: 45 },
    { rating: '4 Stars', count: 32, percentage: 32 },
    { rating: '3 Stars', count: 15, percentage: 15 },
    { rating: '2 Stars', count: 6, percentage: 6 },
    { rating: '1 Star', count: 2, percentage: 2 }
  ],
  categoryRatings: [
    { category: 'Algebra', avgRating: 4.2, count: 68 },
    { category: 'Geometry', avgRating: 4.0, count: 45 },
    { category: 'Calculus', avgRating: 3.8, count: 32 },
    { category: 'Statistics', avgRating: 4.1, count: 23 }
  ],
  recentFeedback: [
    {
      id: '1',
      rating: 5,
      comment: 'The step-by-step guidance was excellent! Really helped me understand the process.',
      timestamp: '2024-01-07 14:30',
      category: 'Algebra',
      problemId: 'p123'
    },
    {
      id: '2',
      rating: 4,
      comment: 'Good hints, but could use more detailed explanations for complex problems.',
      timestamp: '2024-01-07 13:15',
      category: 'Calculus',
      problemId: 'p124'
    },
    {
      id: '3',
      rating: 5,
      comment: 'Perfect balance of guidance without giving away the answer!',
      timestamp: '2024-01-07 12:45',
      category: 'Geometry',
      problemId: 'p125'
    }
  ]
}

const COLORS = ['#22c55e', '#3b82f6', '#f59e0b', '#ef4444', '#6b7280']

export function FeedbackAnalysisPage() {
  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <h1 className="text-3xl font-bold text-gray-900">Feedback Analysis</h1>
        <div className="flex items-center space-x-3">
          <button className="inline-flex items-center px-4 py-2 border border-gray-300 rounded-md text-sm font-medium text-gray-700 bg-white hover:bg-gray-50">
            <Filter className="h-4 w-4 mr-2" />
            Filter
          </button>
          <select className="px-3 py-2 border border-gray-300 rounded-md text-sm focus:outline-none focus:ring-2 focus:ring-primary-500">
            <option>Last 7 days</option>
            <option>Last 30 days</option>
            <option>Last 3 months</option>
          </select>
        </div>
      </div>

      {/* Summary Stats */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-6">
        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <MessageSquare className="h-8 w-8 text-primary-600" />
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Total Feedback</p>
              <p className="text-2xl font-bold text-gray-900">168</p>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <Star className="h-8 w-8 text-yellow-500" />
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Average Rating</p>
              <p className="text-2xl font-bold text-gray-900">4.2</p>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <TrendingUp className="h-8 w-8 text-green-500" />
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Satisfaction Rate</p>
              <p className="text-2xl font-bold text-gray-900">91%</p>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <MessageSquare className="h-8 w-8 text-blue-500" />
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Response Rate</p>
              <p className="text-2xl font-bold text-gray-900">76%</p>
            </div>
          </div>
        </div>
      </div>

      {/* Charts */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <div className="bg-white rounded-lg border p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Rating Distribution</h3>
          <ResponsiveContainer width="100%" height={300}>
            <PieChart>
              <Pie
                data={mockFeedbackData.ratingDistribution}
                cx="50%"
                cy="50%"
                labelLine={false}
                label={({ rating, percentage }) => `${rating}: ${percentage}%`}
                outerRadius={80}
                fill="#8884d8"
                dataKey="count"
              >
                {mockFeedbackData.ratingDistribution.map((entry, index) => (
                  <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                ))}
              </Pie>
              <Tooltip />
            </PieChart>
          </ResponsiveContainer>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Average Rating by Category</h3>
          <ResponsiveContainer width="100%" height={300}>
            <BarChart data={mockFeedbackData.categoryRatings}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="category" />
              <YAxis domain={[0, 5]} />
              <Tooltip />
              <Bar dataKey="avgRating" fill="#3b82f6" />
            </BarChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Recent Feedback */}
      <div className="bg-white rounded-lg border">
        <div className="px-6 py-4 border-b">
          <h3 className="text-lg font-semibold text-gray-900">Recent Feedback</h3>
        </div>
        <div className="divide-y divide-gray-200">
          {mockFeedbackData.recentFeedback.map((feedback) => (
            <div key={feedback.id} className="p-6">
              <div className="flex items-start justify-between">
                <div className="flex-1">
                  <div className="flex items-center space-x-2 mb-2">
                    <div className="flex items-center">
                      {[...Array(5)].map((_, i) => (
                        <Star
                          key={i}
                          className={`h-4 w-4 ${
                            i < feedback.rating
                              ? 'text-yellow-400 fill-current'
                              : 'text-gray-300'
                          }`}
                        />
                      ))}
                    </div>
                    <span className="text-sm text-gray-500">
                      {feedback.category} • {feedback.timestamp}
                    </span>
                  </div>
                  <p className="text-gray-900">{feedback.comment}</p>
                </div>
                <div className="ml-4 text-sm text-gray-500">
                  Problem #{feedback.problemId}
                </div>
              </div>
            </div>
          ))}
        </div>
        <div className="px-6 py-4 border-t bg-gray-50">
          <button className="text-primary-600 hover:text-primary-700 text-sm font-medium">
            View all feedback →
          </button>
        </div>
      </div>

      {/* Insights */}
      <div className="bg-white rounded-lg border p-6">
        <h3 className="text-lg font-semibold text-gray-900 mb-4">Key Insights</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          <div className="space-y-3">
            <h4 className="font-medium text-gray-900">Positive Trends</h4>
            <ul className="space-y-2 text-sm text-gray-600">
              <li className="flex items-start">
                <span className="text-green-500 mr-2">•</span>
                Step-by-step guidance consistently rated highly (4.5+ avg)
              </li>
              <li className="flex items-start">
                <span className="text-green-500 mr-2">•</span>
                91% of users find the hints helpful without being too revealing
              </li>
              <li className="flex items-start">
                <span className="text-green-500 mr-2">•</span>
                Algebra problems show strongest engagement (avg 4.2 rating)
              </li>
            </ul>
          </div>
          <div className="space-y-3">
            <h4 className="font-medium text-gray-900">Areas for Improvement</h4>
            <ul className="space-y-2 text-sm text-gray-600">
              <li className="flex items-start">
                <span className="text-yellow-500 mr-2">•</span>
                Calculus problems need more detailed explanations
              </li>
              <li className="flex items-start">
                <span className="text-yellow-500 mr-2">•</span>
                Some users want more context for advanced concepts
              </li>
              <li className="flex items-start">
                <span className="text-yellow-500 mr-2">•</span>
                Response time could be faster for complex problems
              </li>
            </ul>
          </div>
        </div>
      </div>
    </div>
  )
}