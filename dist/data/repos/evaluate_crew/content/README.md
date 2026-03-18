# Evaluate Crew: Multi-Agent System for Educational Assessment and Improvement

Evaluate Crew is a sophisticated multi-agent system designed to facilitate education through quiz generation, performance analysis, and quality assurance. Built on the CrewAI platform, this system leverages advanced agent configurations to create a sequential process that enhances the learning experience.

## System Overview

- **Quiz Bot (Agent 1)**: Generates a quiz based on a user-selected topic. This agent assesses the user's understanding by crafting questions that test various aspects of the topic, providing scores and performance feedback.
  
- **Analyze Bot (Agent 2)**: Analyzes the quiz results to determine the user's performance and identifies areas needing improvement. It then crafts a detailed learning roadmap from basic to advanced levels, ensuring a structured and efficient learning path.
  
- **Quality Agent (Agent 3)**: Reviews the learning roadmap for accuracy and relevance. If the roadmap meets quality standards, it is finalized and documented in a markdown file; otherwise, it is sent back to the Analyze Bot for revisions.

## Installation

Ensure you have Python installed and follow these steps to set up the Evaluate Crew on your local system:

```bash
git clone https://github.com/yourusername/evaluate-crew.git
cd evaluate-crew
pip install -r requirements.txt
```

### Usage
To run the system, navigate to the project directory and execute the main application:

```bash
python run_app.py
```
The user interface will prompt you to select a topic for the quiz. Once selected, the system proceeds through the agents as follows:

Quiz Bot generates the quiz and collects user responses.
Analyze Bot processes the results and drafts a learning roadmap.
Quality Agent assesses the roadmap and either approves it for documentation or returns it for adjustment.

### Contributing
Contributions to the Evaluate Crew are welcome! If you have suggestions for improvement or new features, please fork the repository, make your changes, and submit a pull request.

### Configuration
Each agent's behavior is defined by YAML configuration files, which specify roles, goals, and operational details. These configurations ensure that the agents perform their tasks effectively and in alignment with the overall objectives of the system.
