// Global error handling middleware for Express
import { Request, Response, NextFunction } from 'express';

export const errorHandler = (
  error: Error,
  req: Request,
  res: Response,
  next: NextFunction
): void => {
  console.error('Error:', error.message);
  console.error('Stack:', error.stack);

  if (res.headersSent) {
    return next(error);
  }

  const statusCode = (error as any).statusCode || 500;
  const message = error.message || 'Internal Server Error';

  res.status(statusCode).json({
    error: message,
    ...(process.env.NODE_ENV === 'development' && { stack: error.stack })
  });
};