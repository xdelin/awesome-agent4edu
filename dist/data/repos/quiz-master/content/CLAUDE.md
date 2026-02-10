# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a vanilla JavaScript quiz application that supports both manual quiz loading and AI-generated quiz creation via backend APIs. The app features user authentication, an interactive quiz interface, and integration with AWS Lambda endpoints for quiz generation and storage.

## Development Commands

Since this is a vanilla HTML/CSS/JavaScript project, there are no build scripts or package managers. Development is done by:

- **Local Development**: Open `index.html` directly in a browser or serve via any local HTTP server
- **Testing**: No automated test framework - manual testing through the browser interface

## Architecture Overview

### Authentication System (`AuthManager` class in script.js:238-343)
- Session-based authentication using `sessionStorage`
- Hardcoded user credentials in the `users` object
- Protects all quiz functionality behind login
- AWS lambda code does its own auth.

### Quiz Management (script.js:1-235)
- **Manual Quiz Loading**: Users can paste JSON quiz data into textarea
- **AI-Generated Quizzes**: Form submission calls backend API to generate questions
- **Quiz State**: Manages current question, score, user answers, and progress

### Navigation & UI
- **Side Menu**: Toggleable navigation with overlay
- **State Management**: Shows/hides different views (login, quiz form, quiz questions)
- **Responsive Design**: CSS handles mobile and desktop layouts

### Backend Integration
- **Configuration**: `config.json` contains AWS Lambda URLs for API endpoints
- **API Calls**: 
  - `QUIZ_FETCH_URL`: Generates quiz questions from prompt
  - `QUIZ_CREATE_URL`: Saves completed quiz data
- **Authentication**: Uses Basic Auth with username/password from session

## Key Files

- `script.js`: Contains all application logic including AuthManager class, quiz functionality, and API integration
- `index.html`: Single page containing login form, side menu, quiz creation form, and quiz interface
- `config.json`: API endpoint configuration for AWS Lambda functions
- `quiz.json`: Sample quiz data used as default content
- `styles.css`: Complete styling for all UI components

## Important Implementation Details

### Quiz Data Structure
```javascript
[
  {
    "question": "Question text",
    "options": ["Option 1", "Option 2", "Option 3", "Option 4"],
    "correct": 0  // Index of correct answer
  }
]
```

### Authentication Flow
1. User logs in via login form
2. Credentials validated against hardcoded user list
3. Session state stored in `sessionStorage`
4. All quiz functions wrapped with auth checks

### API Integration Pattern
- Uses `fetch` with Basic Auth headers
- Error handling with try/catch blocks
- Response data parsing for quiz generation
- Session storage for user credentials across API calls

## Configuration Notes

- Backend API URLs are configured in `config.json`
- Default fallback URLs point to `localhost:3000` if config loading fails
- User credentials are hardcoded in the `AuthManager.users` object