// Routes for math problem solving and step-by-step guidance
import { Router } from 'express';
import { MathController } from '../controllers/MathController';

const router = Router();
const mathController = new MathController();

router.post('/solve', mathController.solveProblem);
router.post('/hint', mathController.getHint);
router.post('/check-step', mathController.checkStep);

export { router as mathRoutes };