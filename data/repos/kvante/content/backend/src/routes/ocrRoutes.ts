// Routes for OCR image processing and text extraction
import { Router } from 'express';
import multer from 'multer';
import { OCRController } from '../controllers/OCRController';

const router = Router();
const ocrController = new OCRController();

const upload = multer({
  storage: multer.memoryStorage(),
  limits: { fileSize: 5 * 1024 * 1024 }, // 5MB limit
  fileFilter: (req, file, cb) => {
    if (file.mimetype.startsWith('image/')) {
      cb(null, true);
    } else {
      cb(new Error('Only image files are allowed'));
    }
  }
});

router.post('/extract', upload.single('image'), ocrController.extractText);

export { router as ocrRoutes };