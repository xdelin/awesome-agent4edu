// Routes for collecting and managing user feedback
import { Router } from 'express';
import { FeedbackController } from '../controllers/FeedbackController';

const router = Router();
const feedbackController = new FeedbackController();

router.post('/', feedbackController.submitFeedback);
router.get('/', feedbackController.getFeedback);

export { router as feedbackRoutes };