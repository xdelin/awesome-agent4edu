# !pip install anthropic docx2txt PyPDF2 spacy nltk
# !python -m spacy download en_core_web_sm

import os
import re
import json
import anthropic
import docx2txt
import spacy
from PyPDF2 import PdfReader
from nltk.corpus import stopwords
from spacy.matcher import Matcher
import nltk
from google.colab import files

nltk.download('stopwords')
nltk.download('punkt')

nlp = spacy.load('en_core_web_sm')
matcher = Matcher(nlp.vocab)

api_key = "your-api-key"

STOPWORDS = set(stopwords.words('english'))
EDUCATION = [
    'BE', 'B.E.', 'B.E', 'BS', 'B.S',
    'ME', 'M.E', 'M.E.', 'M.B.A', 'MBA', 'MS', 'M.S',
    'BTECH', 'B.TECH', 'M.TECH', 'MTECH',
    'SSLC', 'SSC', 'HSC', 'CBSE', 'ICSE', 'X', 'XII'
]
SKILLS = [
    'python', 'java', 'machine learning', 'data analysis', 'sql',
    'html', 'css', 'javascript', 'react', 'node.js', 'aws', 'docker',
    'kubernetes', 'c++', 'c#', 'git', 'linux', 'cloud computing', 'devops'
]

def pdftotext(file_path):
    try:
        reader = PdfReader(file_path)
        text = ""
        for page in reader.pages:
            if page.extract_text():
                text += page.extract_text()
        return text
    except Exception as e:
        print(f"Error reading PDF: {e}")
        return ""

def doctotext(file_path):
    try:
        text = docx2txt.process(file_path)
        return ' '.join([line.strip() for line in text.split('\n') if line.strip()])
    except Exception as e:
        print(f"Error reading DOCX: {e}")
        return ""

def extract_name(resume_text):
    nlp_text = nlp(resume_text)
    pattern = [{'POS': 'PROPN'}, {'POS': 'PROPN'}]
    matcher.add('NAME', [pattern])
    matches = matcher(nlp_text)
    for match_id, start, end in matches:
        span = nlp_text[start:end]
        return span.text
    return None

def extract_education(resume_text):
    nlp_text = nlp(resume_text)
    edu = {}
    for sent in nlp_text.sents:
        for word in sent.text.split():
            word_clean = re.sub(r'[?|$|.|!|,]', r'', word)
            if word_clean.upper() in EDUCATION and word_clean not in STOPWORDS:
                edu[word_clean] = sent.text.strip()
    return list(edu.keys())

def extract_skills(text):
    nlp_text = nlp(text)
    tokens = [token.text.lower() for token in nlp_text if not token.is_stop]
    return list(set(token for token in tokens if token in SKILLS))

def extract_mobile_number(resume_text):
    phone_pattern = r'\b(?:\+?91[-.\s]?)?\d{10}\b|\b\d{3}[-.\s]?\d{3}[-.\s]?\d{4}\b'
    phone = re.findall(phone_pattern, resume_text)
    return phone[0] if phone else None

def extract_email_addresses(resume_text):
    email_pattern = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}\b'
    return re.findall(email_pattern, resume_text)

def extract_keywords(text):
    nlp_text = nlp(text.lower())
    keywords = [token.text for token in nlp_text if token.text in SKILLS]
    return list(set(keywords))

def calculate_match_percentage(jd_keywords, resume_keywords):
    if not jd_keywords:
        return 0
    match_count = len(set(jd_keywords) & set(resume_keywords))
    return (match_count / len(jd_keywords)) * 100

def query_claude(candidates_data, jd_text, query):
    prompt = f"""
    Job Description:
    {jd_text}

    Candidates Data:
    {json.dumps(candidates_data, indent=4)}

    Question: {query}
    Provide a detailed answer based on the above information.
    """
    try:
        client = anthropic.Client(api_key=api_key)
        full_prompt = f"{anthropic.HUMAN_PROMPT}{prompt}{anthropic.AI_PROMPT}"
        response = client.completions.create(
            model="claude-2",
            prompt=full_prompt,
            max_tokens_to_sample=500
        )
        return response.completion.strip()
    except Exception as e:
        return f"Error querying Claude AI: {e}"

def main():
    print("Please upload the Job Description PDF or DOCX file.")
    jd_file = files.upload()
    jd_file_path = next(iter(jd_file))
    if jd_file_path.endswith('.pdf'):
        jd_text = pdftotext(jd_file_path)
    elif jd_file_path.endswith('.docx'):
        jd_text = doctotext(jd_file_path)
    else:
        print("Unsupported JD file type!")
        return
    jd_keywords = extract_keywords(jd_text)
    print("Please upload the resumes (PDF or DOCX files).")
    resumes = files.upload()
    resume_file_paths = list(resumes.keys())
    candidates_data = []
    for resume_file_path in resume_file_paths:
        if resume_file_path.endswith('.pdf'):
            resume_text = pdftotext(resume_file_path)
        elif resume_file_path.endswith('.docx'):
            resume_text = doctotext(resume_file_path)
        else:
            print(f"Unsupported resume file type: {resume_file_path}")
            continue
        resume_keywords = extract_keywords(resume_text)
        match_percentage = calculate_match_percentage(jd_keywords, resume_keywords)
        candidate_data = {
            "Name": extract_name(resume_text),
            "Education": extract_education(resume_text),
            "Skills": extract_skills(resume_text),
            "Mobile Number": extract_mobile_number(resume_text),
            "Email Addresses": extract_email_addresses(resume_text),
            "Match Percentage": match_percentage
        }
        candidates_data.append(candidate_data)
    print("\n--- Extracted Candidate Data ---")
    for candidate in candidates_data:
        print(json.dumps(candidate, indent=4))
    print("\n--- Query System ---")
    while True:
        query = input("Ask a question or type 'exit' to quit: ").strip()
        if query.lower() == "exit":
            break
        answer = query_claude(candidates_data, jd_text, query)
        print("\n--- Answer ---")
        print(answer)

if _name_ == "_main_":
    main()