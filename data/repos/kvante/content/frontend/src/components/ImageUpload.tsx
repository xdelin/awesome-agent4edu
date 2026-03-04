// Component for uploading and processing images with OCR
import { useState, useRef } from 'react'
import { Camera, Upload } from 'lucide-react'
import { useOCR } from '../hooks/useOCR'

interface ImageUploadProps {
  onTextExtracted: (text: string) => void
}

export function ImageUpload({ onTextExtracted }: ImageUploadProps) {
  const [selectedFile, setSelectedFile] = useState<File | null>(null)
  const [previewUrl, setPreviewUrl] = useState<string | null>(null)
  const fileInputRef = useRef<HTMLInputElement>(null)
  const { extractText, isLoading } = useOCR()

  const handleFileSelect = (file: File) => {
    setSelectedFile(file)
    const url = URL.createObjectURL(file)
    setPreviewUrl(url)
  }

  const handleFileInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0]
    if (file) {
      handleFileSelect(file)
    }
  }

  const handleUploadClick = () => {
    fileInputRef.current?.click()
  }

  const handleExtractText = async () => {
    if (selectedFile) {
      const text = await extractText(selectedFile)
      if (text) {
        onTextExtracted(text)
      }
    }
  }

  return (
    <div className="space-y-4">
      <input
        ref={fileInputRef}
        type="file"
        accept="image/*"
        onChange={handleFileInputChange}
        className="hidden"
        capture="environment"
      />
      
      <div className="border-2 border-dashed border-gray-300 rounded-lg p-6 text-center">
        {previewUrl ? (
          <div className="space-y-4">
            <img 
              src={previewUrl} 
              alt="Preview" 
              className="max-w-full h-48 mx-auto object-contain rounded"
            />
            <button
              onClick={handleExtractText}
              disabled={isLoading}
              className="bg-primary-600 text-white px-4 py-2 rounded hover:bg-primary-700 disabled:opacity-50"
            >
              {isLoading ? 'Extracting...' : 'Extract Text'}
            </button>
          </div>
        ) : (
          <div className="space-y-4">
            <Camera className="h-12 w-12 text-gray-400 mx-auto" />
            <div>
              <button
                onClick={handleUploadClick}
                className="bg-primary-600 text-white px-4 py-2 rounded hover:bg-primary-700 inline-flex items-center"
              >
                <Upload className="h-4 w-4 mr-2" />
                Upload Image
              </button>
            </div>
            <p className="text-sm text-gray-500">
              Upload a photo of your math problem
            </p>
          </div>
        )}
      </div>
    </div>
  )
}