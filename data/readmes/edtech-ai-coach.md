# AI-Powered Student Performance Feedback System 📝

## Features
- Parses JSON data, extracting subject- and chapter-level accuracy, time, and difficulty.  
- Generates AI feedback with:  
  - Personalized intro  
  - Performance breakdown  
  - Time vs Accuracy insights (forced to avoid “N/A”)  
  - Actionable suggestions  
- Displays charts in Streamlit UI.  
- Generates downloadable PDF reports with all feedback and visuals.
  
![image](https://github.com/user-attachments/assets/66d4326e-c38f-466e-aa4f-8ab31edef055)
![image](https://github.com/user-attachments/assets/358c9097-58b4-4d5b-a5c4-379775103ae8)
![image](https://github.com/user-attachments/assets/9f20b1e1-6356-482c-9a9f-2f2e9a89c1d4)
![image](https://github.com/user-attachments/assets/6019d171-ae6f-4581-9b25-73033d3c2478)

## Tech Stack
- Python 3.x  
- Streamlit (UI)  
- Matplotlib (charts)  
- Gemini API (LLM)  
- PDF generation libraries (ReportLab/FPDF)  
- JSON parsing

## Prompt Engineering
- Uses structured prompts with examples to ensure human-like, personalized, and motivating feedback.  
- Emphasizes avoiding empty insights by forcing time-vs-accuracy analysis.

## Running the App
1. Install dependencies: `pip install -r requirements.txt`  
2. Set API keys for LLM access in `generate_feedback.py`  
3. Run: `streamlit run app.py`  
4. Upload/select JSON test data and generate feedback.  
5. View AI feedback, charts, and download PDF report.

---


Thank you!
