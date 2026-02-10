// Main App component with routing configuration
import { Routes, Route } from 'react-router-dom'
import { HomePage } from './pages/HomePage'
import { ProblemSolverPage } from './pages/ProblemSolverPage'
import { Layout } from './components/Layout'

function App() {
  return (
    <Layout>
      <Routes>
        <Route path="/" element={<HomePage />} />
        <Route path="/solve" element={<ProblemSolverPage />} />
      </Routes>
    </Layout>
  )
}

export default App