// Main dashboard App component with routing
import { Routes, Route } from 'react-router-dom'
import { Layout } from './components/Layout'
import { OverviewPage } from './pages/OverviewPage'
import { PromptTestingPage } from './pages/PromptTestingPage'
import { FeedbackAnalysisPage } from './pages/FeedbackAnalysisPage'
import { PerformanceMetricsPage } from './pages/PerformanceMetricsPage'

function App() {
  return (
    <Layout>
      <Routes>
        <Route path="/" element={<OverviewPage />} />
        <Route path="/prompt-testing" element={<PromptTestingPage />} />
        <Route path="/feedback" element={<FeedbackAnalysisPage />} />
        <Route path="/metrics" element={<PerformanceMetricsPage />} />
      </Routes>
    </Layout>
  )
}

export default App