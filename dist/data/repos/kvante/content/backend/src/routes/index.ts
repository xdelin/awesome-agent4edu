// Main API routes configuration
import { Router } from 'express';
import { mathRoutes } from './mathRoutes';
import { ocrRoutes } from './ocrRoutes';
import { feedbackRoutes } from './feedbackRoutes';

const router = Router();

router.use('/math', mathRoutes);
router.use('/ocr', ocrRoutes);
router.use('/feedback', feedbackRoutes);

export { router as apiRoutes };