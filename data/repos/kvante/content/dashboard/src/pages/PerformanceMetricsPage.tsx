// Performance metrics and system monitoring page
import { LineChart, Line, AreaChart, Area, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts'
import { Clock, Zap, AlertTriangle, CheckCircle2 } from 'lucide-react'

const mockMetrics = {
  responseTime: [
    { time: '00:00', avg: 1.2, p95: 2.1, p99: 3.8 },
    { time: '04:00', avg: 0.9, p95: 1.8, p99: 3.2 },
    { time: '08:00', avg: 1.5, p95: 2.8, p99: 4.5 },
    { time: '12:00', avg: 2.1, p95: 3.5, p99: 5.2 },
    { time: '16:00', avg: 1.8, p95: 3.0, p99: 4.8 },
    { time: '20:00', avg: 1.3, p95: 2.3, p99: 3.9 },
  ],
  systemHealth: [
    { time: '00:00', cpu: 45, memory: 62, requests: 120 },
    { time: '04:00', cpu: 32, memory: 58, requests: 85 },
    { time: '08:00', cpu: 68, memory: 71, requests: 180 },
    { time: '12:00', cpu: 82, memory: 78, requests: 250 },
    { time: '16:00', cpu: 75, memory: 74, requests: 220 },
    { time: '20:00', cpu: 55, memory: 65, requests: 160 },
  ],
  errorRates: [
    { time: '00:00', rate: 0.2 },
    { time: '04:00', rate: 0.1 },
    { time: '08:00', rate: 0.5 },
    { time: '12:00', rate: 1.2 },
    { time: '16:00', rate: 0.8 },
    { time: '20:00', rate: 0.3 },
  ]
}

export function PerformanceMetricsPage() {
  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <h1 className="text-3xl font-bold text-gray-900">Performance Metrics</h1>
        <div className="flex items-center space-x-3">
          <div className="flex items-center space-x-2">
            <div className="h-3 w-3 bg-green-500 rounded-full"></div>
            <span className="text-sm text-gray-600">System Healthy</span>
          </div>
          <select className="px-3 py-2 border border-gray-300 rounded-md text-sm focus:outline-none focus:ring-2 focus:ring-primary-500">
            <option>Last 24 hours</option>
            <option>Last 7 days</option>
            <option>Last 30 days</option>
          </select>
        </div>
      </div>

      {/* Key Metrics */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-6">
        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <Clock className="h-8 w-8 text-blue-600" />
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Avg Response Time</p>
              <p className="text-2xl font-bold text-gray-900">1.4s</p>
              <p className="text-sm text-green-600">-12% vs yesterday</p>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <Zap className="h-8 w-8 text-yellow-600" />
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Throughput</p>
              <p className="text-2xl font-bold text-gray-900">156 req/min</p>
              <p className="text-sm text-green-600">+8% vs yesterday</p>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <AlertTriangle className="h-8 w-8 text-red-600" />
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Error Rate</p>
              <p className="text-2xl font-bold text-gray-900">0.3%</p>
              <p className="text-sm text-red-600">+0.1% vs yesterday</p>
            </div>
          </div>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <div className="flex items-center">
            <CheckCircle2 className="h-8 w-8 text-green-600" />
            <div className="ml-4">
              <p className="text-sm font-medium text-gray-500">Uptime</p>
              <p className="text-2xl font-bold text-gray-900">99.8%</p>
              <p className="text-sm text-green-600">SLA: 99.5%</p>
            </div>
          </div>
        </div>
      </div>

      {/* Response Time Chart */}
      <div className="bg-white rounded-lg border p-6">
        <h3 className="text-lg font-semibold text-gray-900 mb-4">Response Time Distribution</h3>
        <ResponsiveContainer width="100%" height={300}>
          <LineChart data={mockMetrics.responseTime}>
            <CartesianGrid strokeDasharray="3 3" />
            <XAxis dataKey="time" />
            <YAxis />
            <Tooltip />
            <Line type="monotone" dataKey="avg" stroke="#3b82f6" strokeWidth={2} name="Average" />
            <Line type="monotone" dataKey="p95" stroke="#f59e0b" strokeWidth={2} name="95th Percentile" />
            <Line type="monotone" dataKey="p99" stroke="#ef4444" strokeWidth={2} name="99th Percentile" />
          </LineChart>
        </ResponsiveContainer>
      </div>

      {/* System Health */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <div className="bg-white rounded-lg border p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">System Resources</h3>
          <ResponsiveContainer width="100%" height={300}>
            <AreaChart data={mockMetrics.systemHealth}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="time" />
              <YAxis />
              <Tooltip />
              <Area type="monotone" dataKey="cpu" stackId="1" stroke="#3b82f6" fill="#3b82f6" fillOpacity={0.6} name="CPU %" />
              <Area type="monotone" dataKey="memory" stackId="2" stroke="#22c55e" fill="#22c55e" fillOpacity={0.6} name="Memory %" />
            </AreaChart>
          </ResponsiveContainer>
        </div>

        <div className="bg-white rounded-lg border p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Error Rate</h3>
          <ResponsiveContainer width="100%" height={300}>
            <AreaChart data={mockMetrics.errorRates}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="time" />
              <YAxis />
              <Tooltip />
              <Area type="monotone" dataKey="rate" stroke="#ef4444" fill="#ef4444" fillOpacity={0.6} name="Error Rate %" />
            </AreaChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Alerts & Issues */}
      <div className="bg-white rounded-lg border">
        <div className="px-6 py-4 border-b">
          <h3 className="text-lg font-semibold text-gray-900">Recent Alerts & Issues</h3>
        </div>
        <div className="divide-y divide-gray-200">
          <div className="p-6 flex items-start space-x-4">
            <div className="flex-shrink-0">
              <div className="h-2 w-2 bg-yellow-500 rounded-full mt-2"></div>
            </div>
            <div className="flex-1">
              <p className="text-sm font-medium text-gray-900">High CPU usage detected</p>
              <p className="text-sm text-gray-500">CPU usage spiked to 85% at 12:15 PM</p>
              <p className="text-xs text-gray-400">2 hours ago</p>
            </div>
            <div className="text-sm text-yellow-600">Warning</div>
          </div>
          
          <div className="p-6 flex items-start space-x-4">
            <div className="flex-shrink-0">
              <div className="h-2 w-2 bg-green-500 rounded-full mt-2"></div>
            </div>
            <div className="flex-1">
              <p className="text-sm font-medium text-gray-900">Database connection restored</p>
              <p className="text-sm text-gray-500">Connection timeout issues resolved</p>
              <p className="text-xs text-gray-400">4 hours ago</p>
            </div>
            <div className="text-sm text-green-600">Resolved</div>
          </div>
          
          <div className="p-6 flex items-start space-x-4">
            <div className="flex-shrink-0">
              <div className="h-2 w-2 bg-blue-500 rounded-full mt-2"></div>
            </div>
            <div className="flex-1">
              <p className="text-sm font-medium text-gray-900">Deployment completed</p>
              <p className="text-sm text-gray-500">New prompt templates deployed successfully</p>
              <p className="text-xs text-gray-400">6 hours ago</p>
            </div>
            <div className="text-sm text-blue-600">Info</div>
          </div>
        </div>
      </div>

      {/* Performance Summary */}
      <div className="bg-white rounded-lg border p-6">
        <h3 className="text-lg font-semibold text-gray-900 mb-4">Performance Summary</h3>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <div className="text-center">
            <div className="text-3xl font-bold text-green-600">98.2%</div>
            <div className="text-sm text-gray-500">Requests served successfully</div>
          </div>
          <div className="text-center">
            <div className="text-3xl font-bold text-blue-600">1.4s</div>
            <div className="text-sm text-gray-500">Average response time</div>
          </div>
          <div className="text-center">
            <div className="text-3xl font-bold text-purple-600">2.1M</div>
            <div className="text-sm text-gray-500">Total requests processed</div>
          </div>
        </div>
        <div className="mt-6 p-4 bg-green-50 rounded-lg">
          <div className="flex">
            <CheckCircle2 className="h-5 w-5 text-green-500 flex-shrink-0 mt-0.5" />
            <div className="ml-3">
              <p className="text-sm font-medium text-green-800">System Performance Excellent</p>
              <p className="text-sm text-green-700">All metrics within acceptable ranges. No action required.</p>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}