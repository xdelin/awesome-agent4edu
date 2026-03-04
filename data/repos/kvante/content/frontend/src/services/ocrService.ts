// Service for OCR text extraction API calls
import axios from 'axios'

const api = axios.create({
  baseURL: '/api/ocr',
  timeout: 30000,
})

export const ocrService = {
  async extractText(file: File): Promise<string> {
    const formData = new FormData()
    formData.append('image', file)
    
    const response = await api.post('/extract', formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
    })
    
    return response.data.text
  }
}