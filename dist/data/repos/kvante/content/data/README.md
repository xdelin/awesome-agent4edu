# Data Directory

This directory contains local JSON storage for the Kvante application.

## Structure

```
data/
├── feedback/          # User feedback and ratings
├── sessions/          # User session data
├── problems/          # Math problem data and solutions
└── logs/             # Application logs
```

## Data Files

All data files are stored as JSON and are gitignored for privacy and size reasons. The directory structure is preserved through `.gitkeep` files.

### Feedback Data
- User ratings (1-5 stars)
- Comments and suggestions
- Session and problem associations
- Timestamps

### Session Data
- User session tracking
- Problem progression
- Performance metrics
- Learning analytics

### Problem Data
- Math problem definitions
- AI-generated solutions
- Step-by-step breakdowns
- Difficulty classifications

### Logs
- Application performance logs
- Error tracking
- API request/response logs
- User interaction logs

## Usage

The backend services will automatically create and manage JSON files in these directories as users interact with the application. No manual intervention is typically required.