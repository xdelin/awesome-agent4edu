from typing import Dict, Any

def generate_feedback_prompt(summary: Dict[str, Any]) -> str:
    """
    Generates a highly personalized, motivating, and constructive feedback prompt
    for an LLM, robust to missing data, and effective for all performance levels.
    Includes a strict instruction to write a meaningful Time vs Accuracy insight
    based on raw data, avoiding 'N/A' or skipping.
    """

    overall_acc   = summary.get("overall_accuracy", 0)
    total_score   = summary.get("total_score", 0)
    max_score     = summary.get("max_score", "N/A")
    subjects      = summary.get("subjects", [])
    chapters      = summary.get("chapters", [])
    suggestions   = summary.get("actionable_suggestions", [])
    student_name  = summary.get("student_name", "the student")

    taraw = summary.get("time_accuracy_raw", {})
    avg_time = taraw.get("avg_time_per_question", None)
    avg_acc = taraw.get("avg_accuracy_percent", None)
    sample_times = taraw.get("sample_times_sec", [])
    sample_correct = taraw.get("sample_correctness", [])

    subj_lines = []
    for s in subjects:
        subj_lines.append(
            f"- **{s.get('name', 'Unknown')}**: "
            f"{s.get('accuracy', 0):.1f}% accuracy, "
            f"{s.get('time_taken', 0)} s spent, "
            f"score {s.get('score', 0)}"
        )
    subj_block = "\n".join(subj_lines) if subj_lines else "_No subject data available._"

    chapters_sorted = sorted(chapters, key=lambda c: c.get("accuracy", 0))
    weakest = chapters_sorted[:3] if chapters_sorted else []
    strongest = chapters_sorted[-3:] if len(chapters_sorted) >= 3 else chapters_sorted

    def _fmt_chapter(c: Dict[str, Any]) -> str:
        return (f"* {c.get('name', 'Chap?')} "
                f"({c.get('difficulty', 'N/A')}): "
                f"{c.get('accuracy', 0):.1f}% acc, "
                f"{c.get('time_taken', 0)} s")

    weak_block = "\n".join(_fmt_chapter(c) for c in weakest) or "_N/A_"
    strong_block = "\n".join(_fmt_chapter(c) for c in strongest) or "_N/A_"

    sugg_block = "\n".join(f"- {s}" for s in suggestions) or "_No suggestions_"

    times_str = ", ".join(str(t) for t in sample_times) if sample_times else "No data"
    correct_str = ", ".join(str(c) for c in sample_correct) if sample_correct else "No data"

    intro_guidance = (
        "Begin your feedback with a warm, empathetic, and encouraging message that "
        "matches the student's overall performance. If the performance is strong, "
        "celebrate achievements. If it is average, encourage and highlight progress. "
        "If it is weak, be supportive, avoid criticism, and focus on growth and potential. "
        "Always sound human, not generic or robotic.\n\n"
        "Provide a clear, insightful, and data-driven analysis of how the student's pacing (time spent per question) "
        "relates to their accuracy, using the sample question times and correctness. "
        "This Time vs Accuracy insight is mandatory and must be meaningful, never left as 'N/A' or omitted.\n\n"
        "Give a detailed breakdown of performance by subjects and chapters, emphasizing areas of strength and improvement.\n\n"
        "Conclude with 2–3 concrete, actionable suggestions tailored to the student's performance profile.\n\n"
        "Use markdown formatting with short paragraphs, bullet points, and tables if needed for clarity."
    )

    prompt_md = f"""
You are an expert, empathetic tutor AI. Write **personalized, motivating, and constructive feedback** for {student_name} based on the performance data below.

{intro_guidance}

---

## Overall Performance
- **Accuracy:** {overall_acc:.1f}%
- **Score:** {total_score} / {max_score}

## Subject Performance
{subj_block}

## Strongest Chapters
{strong_block}

## Weakest Chapters
{weak_block}

## Time vs Accuracy Raw Data
- Average time per question: {avg_time if avg_time is not None else 'No data'} seconds
- Average correctness: {avg_acc if avg_acc is not None else 'No data'}%
- Sample question times (sec): [{times_str}]
- Sample question correctness (1=correct,0=incorrect): [{correct_str}]

---

**Important:** Provide a detailed and insightful analysis of how the student's question response times correlate with their accuracy. Discuss whether they perform better with more or less time, identify any pacing patterns, and suggest how they might optimize their time management based on this data.

## Actionable Suggestions
{sugg_block}
"""
    return prompt_md.strip()
