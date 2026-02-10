import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import os
from generate_feedback import (
    load_single_report,
    process_reports_to_summary,
    generate_feedback_prompt,
    call_gemini_api
)
from fpdf import FPDF
from fpdf.enums import XPos, YPos  
import re
import io

st.set_page_config(page_title="Student Performance Feedback", layout="wide")

# ---------------- Chart helpers ---------------- #
def plot_subject_accuracy(subjects):
    if not subjects:
        st.warning("No subject data to plot.")
        return None
    df = pd.DataFrame(subjects)
    if df.empty or "name" not in df.columns or "accuracy" not in df.columns:
        st.warning("Subject data is missing required fields.")
        return None
    fig, ax = plt.subplots()
    ax.bar(df["name"], df["accuracy"], color="skyblue")
    ax.set_ylim(0, 100)
    ax.set_ylabel("Accuracy (%)")
    ax.set_title("Subject-wise Accuracy")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    return fig

def plot_chapter_time(chapters):
    if not chapters:
        st.warning("No chapter data to plot.")
        return None
    df = pd.DataFrame(chapters)
    if df.empty or "name" not in df.columns or "time_taken" not in df.columns:
        st.warning("Chapter data is missing required fields.")
        return None
    fig, ax = plt.subplots()
    ax.barh(df["name"], df["time_taken"], color="coral")
    ax.set_xlabel("Time Taken (s)")
    ax.set_title("Chapter-wise Time Spent")
    plt.tight_layout()
    return fig

# ---------------- PDF helper ----------------
from fpdf import FPDF
from fpdf.enums import XPos, YPos
import io
import re

class PDF(FPDF):
    def header(self):
        self.set_font("Helvetica", "B", 16)
        self.cell(0, 10, "Student Performance Feedback Report", align="C", new_x=XPos.LMARGIN, new_y=YPos.NEXT)

def sanitize_text(text: str) -> str:
    if not isinstance(text, str):
        text = str(text)
    replacements = {
        "’": "'", "‘": "'", "“": '"', "”": '"', "–": "-", "—": "-",
        "…": "...", "•": "-", "\u2013": "-", "\u2014": "-",
        "é": "e", "á": "a", "à": "a", "ç": "c", "ó": "o", "ö": "o", "ü": "u"
    }
    for k, v in replacements.items():
        text = text.replace(k, v)
    text = re.sub(r"[\x00-\x08\x0B\x0C\x0E-\x1F\x7F]", "", text)
    return text

def parse_feedback_to_pdf(pdf, text):
    text = sanitize_text(text)
    lines = text.split('\n')
    for line in lines:
        line = line.strip()
        if not line:
            continue  
        if line.startswith(("-", "*")):
            pdf.set_font("Helvetica", "", 11)
            pdf.multi_cell(0, 6, "  " + line, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        else:
            pdf.set_font("Helvetica", "", 11)
            pdf.multi_cell(0, 6, line, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

def generate_pdf(summary, feedback_text, figs=None) -> bytes:
    pdf = PDF()
    pdf.add_page()
    pdf.set_auto_page_break(auto=True, margin=15)
    line_height = 8 
    page_width = pdf.w - 2 * pdf.l_margin

    # Overall Accuracy
    overall_acc = summary.get("overall_accuracy", 0)
    pdf.set_font("Helvetica", "", 11)
    pdf.multi_cell(0, line_height, f"Overall Accuracy: {overall_acc:.2f}%", new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    # Subject Table
    col_widths = [page_width * 0.45, page_width * 0.25, page_width * 0.25]
    headers = ["Subject", "Accuracy (%)", "Time (s)"]

    pdf.set_font("Helvetica", "B", 11)
    for i, header in enumerate(headers):
        pdf.cell(col_widths[i], line_height, header, border=1, align="C")
    pdf.ln(line_height)

    pdf.set_font("Helvetica", "", 11)
    subjects = summary.get("subjects", [])
    for s in subjects:
        name = sanitize_text(s.get("name", "Unknown"))[:25]  
        accuracy = s.get("accuracy", 0)
        time_taken = s.get("time_taken", 0)
        pdf.cell(col_widths[0], line_height, name, border=1)
        pdf.cell(col_widths[1], line_height, f"{accuracy:.1f}", border=1, align="R")
        pdf.cell(col_widths[2], line_height, str(time_taken), border=1, align="R", new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    
    # Feedback section
    pdf.ln(1)
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, line_height, "AI-Generated Feedback:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    parse_feedback_to_pdf(pdf, feedback_text)

    # Actionable Suggestions
    pdf.ln(1)
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, line_height, "Actionable Suggestions:", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    suggestions = summary.get("actionable_suggestions", ["No suggestions available"])
    pdf.set_font("Helvetica", "", 11)
    for s in suggestions:
        pdf.multi_cell(0, line_height, f"- {sanitize_text(s)}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    # Charts (figures)
    if figs:
        for i, fig in enumerate(figs):
            if fig is None:
                continue
            pdf.add_page()
            img_buffer = io.BytesIO()
            fig.savefig(img_buffer, format='PNG', bbox_inches='tight')
            img_buffer.seek(0)
           
            img_width, img_height = fig.get_size_inches()
            aspect_ratio = img_height / img_width
            display_width = page_width
            display_height = display_width * aspect_ratio
            pdf.image(img_buffer, x=pdf.l_margin, y=pdf.t_margin, w=display_width, h=display_height)
            img_buffer.close()

    return bytes(pdf.output(dest="S"))


# ---------------- Streamlit main ---------------- #
def main():
    st.title("📊 Student Performance Feedback Generator")

    data_folder = st.text_input("Data folder path:", value="data")

    if not os.path.isdir(data_folder):
        st.error("Invalid folder path. Please enter a valid data folder path.")
        return

    files = [f for f in os.listdir(data_folder) if f.lower().endswith(".json")]
    if not files:
        st.error("No JSON files found in the folder.")
        return

    selected_file = st.selectbox("Select a JSON report file to process:", files)

    if st.button("Generate Feedback"):
        with st.spinner(f"Processing {selected_file} ..."):
            filepath = os.path.join(data_folder, selected_file)
            try:
                report = load_single_report(filepath)
                summary = process_reports_to_summary([report])
                prompt = generate_feedback_prompt(summary)
                feedback_text = call_gemini_api(prompt)
            except Exception as e:
                st.error(f"Error processing file or generating feedback: {e}")
                return

        st.header("✨ Personalized AI Feedback")
        st.markdown(feedback_text)

        st.header("📈 Performance Breakdown")

        fig_sub = plot_subject_accuracy(summary.get("subjects", []))
        if fig_sub:
            st.pyplot(fig_sub)

        fig_ch = plot_chapter_time(summary.get("chapters", []))
        if fig_ch:
            st.pyplot(fig_ch)

        

        st.subheader("💡 Actionable Suggestions")
        for s in summary.get("actionable_suggestions", []):
            st.write(f"- {s}")

        try:
            figs = []
            if fig_sub:
                figs.append(fig_sub)
            if fig_ch:
                figs.append(fig_ch)

            pdf_bytes = generate_pdf(summary, feedback_text, figs)
            st.download_button(
                "📥 Download PDF report",
                data=pdf_bytes,
                file_name=f"student_feedback_{selected_file.replace('.json', '')}.pdf",
                mime="application/pdf",
            )
        except Exception as e:
            st.error(f"Error generating PDF: {e}")

if __name__ == "__main__":
    main()
