// Controller for handling OCR image processing
import { Request, Response } from 'express';
import { OCRService } from '../services/OCRService';

export class OCRController {
  private ocrService: OCRService;

  constructor() {
    this.ocrService = new OCRService();
  }

  extractText = async (req: Request, res: Response): Promise<void> => {
    try {
      if (!req.file) {
        res.status(400).json({ error: 'No image file provided' });
        return;
      }

      const extractedText = await this.ocrService.extractText(req.file.buffer);
      res.json({ text: extractedText });
    } catch (error) {
      console.error('Error extracting text from image:', error);
      res.status(500).json({ error: 'Failed to extract text from image' });
    }
  };
}