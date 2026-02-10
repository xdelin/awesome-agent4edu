// Dashboard overview page with key metrics and summaries
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, LineChart, Line } from 'recharts'
import { TrendingUp, Users, MessageCircle, CheckCircle } from 'lucide-react'

const mockData = {
  dailyUsage: [
    { date: '2024-01-01', sessions: 45, problems: 120 },
    { date: '2024-01-02', sessions: 52, problems: 145 },
    { date: '2024-01-03', sessions: 38, problems: 98 },
    { date: '2024-01-04', sessions: 61, problems: 167 },
    { date: '2024-01-05', sessions: 48, problems: 134 },
    { date: '2024-01-06', sessions: 55, problems: 156 },
    { date: '2024-01-07', sessions: 43, problems: 115 },
  ],
  promptPerformance: [
    { prompt: 'Step-by-Step v1', rating: 4.2, usage: 45 },
    { prompt: 'Step-by-Step v2', rating: 4.6, usage: 32 },
    { prompt: 'Hint Generator v1', rating: 3.8, usage: 28 },
    { prompt: 'Step Checker v1', rating: 4.1, usage: 35 },
  ]
}

export function OverviewPage() {
  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <h1 className="text-3xl font-bold text-gray-900">Dashboard Overview</h1>
        <div className="text-sm text-gray-500">
          Last updated: {new Date().toLocaleString()}
        </div>
      </div>

      {/* Key Metrics */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <div className="flex-shrink-0">
              <Users className="h-8 w-8 text-primary-600" />
            </div>
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Total Sessions</p>
              <p className="text-2xl font-bold text-gray-900">342</p>
              <p className="text-sm text-green-600 flex items-center">
                <TrendingUp className="h-4 w-4 mr-1" />
                +12.5%
              </p>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <div className="flex-shrink-0">
              <CheckCircle className="h-8 w-8 text-success-600" />
            </div>
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Problems Solved</p>
              <p className="text-2xl font-bold text-gray-900">1,235</p>
              <p className="text-sm text-green-600 flex items-center">
                <TrendingUp className="h-4 w-4 mr-1" />
                +8.3%
              </p>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <div className="flex-shrink-0">
              <MessageCircle className="h-8 w-8 text-warning-600" />
            </div>
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Avg. Rating</p>
              <p className="text-2xl font-bold text-gray-900">4.2</p>
              <p className="text-sm text-green-600 flex items-center">
                <TrendingUp className="h-4 w-4 mr-1" />
                +0.3
              </p>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <div className="flex-shrink-0">
              <TrendingUp className="h-8 w-8 text-error-600" />
            </div>
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Success Rate</p>
              <p className="text-2xl font-bold text-gray-900">87%</p>
              <p className="text-sm text-green-600 flex items-center">
                <TrendingUp className="h-4 w-4 mr-1" />
                +5.2%
              </p>
            </div>
          </div>
        </div>
      </div>

      {/* Charts */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <div className="bg-white rounded-lg border p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Daily Usage</h3>
          <ResponsiveContainer width="100%" height={300}>
            <LineChart data={mockData.dailyUsage}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="date" />
              <YAxis />
              <Tooltip />
              <Line type="monotone" dataKey="sessions" stroke="#3b82f6" strokeWidth={2} />
              <Line type="monotone" dataKey="problems" stroke="#22c55e" strokeWidth={2} />
            </LineChart>
          </ResponsiveContainer>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Prompt Performance</h3>
          <ResponsiveContainer width="100%" height={300}>
            <BarChart data={mockData.promptPerformance}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="prompt" />
              <YAxis />
              <Tooltip />
              <Bar dataKey="rating" fill="#3b82f6" />
            </BarChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Recent Activity */}
      <div className="bg-white rounded-lg border">
        <div className="px-6 py-4 border-b">
          <h3 className="text-lg font-semibold text-gray-900">Recent Activity</h3>
        </div>
        <div className="p-6">
          <div className="space-y-4">
            <div className="flex items-center space-x-4">
              <div className="flex-shrink-0">
                <CheckCircle className="h-5 w-5 text-green-500" />
              </div>
              <div className="flex-1">
                <p className="text-sm font-medium text-gray-900">New prompt version deployed</p>
                <p className="text-sm text-gray-500">Step-by-Step Solver v2.1 - 2 hours ago</p>
              </div>
            </div>
            <div className="flex items-center space-x-4">
              <div className="flex-shrink-0">
                <MessageCircle className="h-5 w-5 text-blue-500" />
              </div>
              <div className="flex-1">
                <p className="text-sm font-medium text-gray-900">High satisfaction feedback received</p>
                <p className="text-sm text-gray-500">Average rating: 4.8/5 - 4 hours ago</p>
              </div>
            </div>
            <div className="flex items-center space-x-4">
              <div className="flex-shrink-0">
                <TrendingUp className="h-5 w-5 text-yellow-500" />
              </div>
              <div className="flex-1">
                <p className="text-sm font-medium text-gray-900">Usage spike detected</p>
                <p className="text-sm text-gray-500">150% increase in algebra problems - 6 hours ago</p>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}