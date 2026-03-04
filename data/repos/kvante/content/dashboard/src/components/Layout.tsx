// Dashboard layout with navigation sidebar
import { ReactNode } from 'react'
import { Link, useLocation } from 'react-router-dom'
import { BarChart3, TestTube, MessageSquare, Activity, Home } from 'lucide-react'

interface LayoutProps {
  children: ReactNode
}

const navigation = [
  { name: 'Overview', href: '/', icon: Home },
  { name: 'Prompt Testing', href: '/prompt-testing', icon: TestTube },
  { name: 'Feedback Analysis', href: '/feedback', icon: MessageSquare },
  { name: 'Performance Metrics', href: '/metrics', icon: Activity },
]

export function Layout({ children }: LayoutProps) {
  const location = useLocation()

  return (
    <div className="min-h-screen flex bg-gray-50">
      <div className="flex flex-col w-64 bg-white shadow-sm">
        <div className="flex items-center h-16 px-6 border-b">
          <BarChart3 className="h-8 w-8 text-primary-600" />
          <span className="ml-2 text-xl font-bold text-gray-900">Kvante Dashboard</span>
        </div>
        
        <nav className="flex-1 px-4 py-6 space-y-2">
          {navigation.map((item) => {
            const Icon = item.icon
            const isActive = location.pathname === item.href
            
            return (
              <Link
                key={item.name}
                to={item.href}
                className={`flex items-center px-3 py-2 text-sm font-medium rounded-md transition-colors ${
                  isActive
                    ? 'bg-primary-50 text-primary-700 border-r-2 border-primary-700'
                    : 'text-gray-600 hover:bg-gray-50 hover:text-gray-900'
                }`}
              >
                <Icon className="h-5 w-5 mr-3" />
                {item.name}
              </Link>
            )
          })}
        </nav>
      </div>
      
      <div className="flex-1 flex flex-col overflow-hidden">
        <main className="flex-1 overflow-auto p-6">
          {children}
        </main>
      </div>
    </div>
  )
}