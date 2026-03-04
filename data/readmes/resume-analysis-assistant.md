# Resume Parser and Matcher Using NLP

A Python-based tool for extracting, analyzing, and matching resumes with job descriptions using advanced Natural Language Processing (NLP) techniques. This tool supports PDF and DOCX file formats and provides insights into resumes by extracting key information such as names, education, skills, contact information, and calculating match percentages with job descriptions. It also integrates with Anthropic's Claude AI for answering advanced queries.

---

## Features

- **Resume Parsing**:
  - Extract text from resumes in PDF or DOCX formats.
  - Identify names, education, skills, phone numbers, and email addresses.
- **Job Description Analysis**:
  - Extracts keywords and analyzes job descriptions.
  - Matches skills with resumes to calculate compatibility percentages.
- **AI Integration**:
  - Uses Anthropic's Claude AI to answer detailed queries about candidates and job descriptions.
- **Keyword Matching**:
  - Leverages NLP to extract domain-specific keywords for analysis.

---

## Prerequisites

- Python 3.7 or higher
- Required Python libraries:
  - `anthropic`
  - `docx2txt`
  - `PyPDF2`
  - `spacy`
  - `nltk`

### Installation Steps

1. Clone the repository:
   
git clone https://github.com/<your-username>/resume-parser.git

Navigate to the directory:
   
cd resume-parser

Install dependencies:
pip install -r requirements.txt

Download spaCy's English language model:
python -m spacy download en_core_web_sm

### Usage

Run the Script:
python resume_parser.py

Upload Files:

When prompted, upload the job description in PDF or DOCX format.
Upload one or more resumes in PDF or DOCX formats.

View Results:
Extracted candidate details and match percentages are displayed.
Use the integrated query system to ask questions about the data.

## Key Functionalities
Extract Candidate Details:
Name
Education
Skills
Contact Information (Phone and Email)
Skill Matching:
Compares resume skills with job description keywords.
Calculates a compatibility percentage.

### Query System:

Allows users to interact with data via Anthropic's Claude AI.
Example: "Who is the best candidate for this role?"
Example Output
Candidate Data
{
    "Name": "John Doe",
    "Education": ["B.Tech", "MBA"],
    "Skills": ["python", "data analysis", "sql"],
    "Mobile Number": "9876543210",
    "Email Addresses": ["johndoe@example.com"],
    "Match Percentage": 85.0
}

### Query Example

Question: "Who is the best candidate for this job?"
Answer: "John Doe with a match percentage of 85% is the most suitable candidate."

### Important Notes
API Key:
Replace api_key in the code with your actual Anthropic API key.
Do not upload your API key to public repositories.

### Customization:
Update the EDUCATION and SKILLS lists in the script for specific industries.
