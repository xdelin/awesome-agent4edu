// Utility functions for OCR input validation
export function validateImageBuffer(buffer: Buffer): boolean {
  if (!buffer || buffer.length === 0) {
    return false
  }

  const imageSignatures = [
    { signature: [0xFF, 0xD8, 0xFF], format: 'JPEG' },
    { signature: [0x89, 0x50, 0x4E, 0x47], format: 'PNG' },
    { signature: [0x47, 0x49, 0x46], format: 'GIF' },
    { signature: [0x42, 0x4D], format: 'BMP' },
    { signature: [0x52, 0x49, 0x46, 0x46], format: 'WEBP' }
  ]

  for (const { signature } of imageSignatures) {
    if (signature.every((byte, index) => buffer[index] === byte)) {
      return true
    }
  }

  return false
}

export function validateImageSize(buffer: Buffer, maxSizeMB: number = 5): boolean {
  const maxSizeBytes = maxSizeMB * 1024 * 1024
  return buffer.length <= maxSizeBytes
}

export function getImageFormat(buffer: Buffer): string | null {
  const signatures = [
    { signature: [0xFF, 0xD8, 0xFF], format: 'JPEG' },
    { signature: [0x89, 0x50, 0x4E, 0x47], format: 'PNG' },
    { signature: [0x47, 0x49, 0x46], format: 'GIF' },
    { signature: [0x42, 0x4D], format: 'BMP' },
    { signature: [0x52, 0x49, 0x46, 0x46], format: 'WEBP' }
  ]

  for (const { signature, format } of signatures) {
    if (signature.every((byte, index) => buffer[index] === byte)) {
      return format
    }
  }

  return null
}