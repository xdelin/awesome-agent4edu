# analytics.py
"""
Step 1 – Parse MathonGo test-report JSON and build an 'LLM context' object.

Usage
-----
python -m analytics <path_to_report.json>

Dependencies
------------
pip install pandas bs4
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from bs4 import BeautifulSoup


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _parse_syllabus(html: str) -> Dict[str, List[str]]:
    """Return {subject -> [chapter, …]} from the embedded syllabus HTML."""
    soup = BeautifulSoup(html, "html.parser")
    out: Dict[str, List[str]] = {}
    for h2 in soup.find_all("h2"):
        subject = h2.text.strip()
        chapters = [
            li.text.strip() for li in h2.find_next_sibling("ul").find_all("li")
        ]
        out[subject] = chapters
    return out


def _subject_from_section(title: str) -> str:
    """
    Section titles start with the subject name, e.g. 'Physics Single Correct'.
    """
    return title.split()[0]  # 'Physics', 'Chemistry', 'Maths', ...


# ---------------------------------------------------------------------------
# core extraction
# ---------------------------------------------------------------------------

def build_dataframes(raw: Dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns
    -------
    df_subjects : one row per subject (already present in JSON).
    df_questions : one row per question with subject, chapter, difficulty, etc.
    """
    # --- subject-level ------------------------------------------------------
    df_subjects = pd.json_normalize(raw["subjects"]).rename(
        columns={
            "totalTimeTaken": "time_taken",
            "totalMarkScored": "marks_scored",
            "totalAttempted": "attempted",
            "totalCorrect": "correct",
        }
    )
    # Add readable subject name by parsing syllabus or section titles
    syllabus_map = _parse_syllabus(raw["test"]["syllabus"])
    oid_to_name = {v: k for k, v in syllabus_map.items() for v in []}  # empty – fallback
    df_subjects["subject"] = df_subjects["subjectId.$oid"].map(
        oid_to_name
    )  # may end up NaN

    # --- question-level -----------------------------------------------------
    records = []
    for section in raw["sections"]:
        subject = _subject_from_section(section["sectionId"]["title"])
        for q in section["questions"]:
            qid = q["questionId"]
            records.append(
                {
                    "subject": subject,
                    "chapter": qid["chapters"][0]["title"],
                    "topic": qid["topics"][0]["title"],
                    "concept": qid["concepts"][0]["title"],
                    "difficulty": qid.get("level", "unknown"),
                    "time_taken": q["timeTaken"],
                    "correct": q["markedOptions"][0]["isCorrect"]
                    if q["markedOptions"]
                    else False,
                }
            )
    df_questions = pd.DataFrame.from_records(records)
    return df_subjects, df_questions


def summarize_subjects(df_subjects: pd.DataFrame) -> List[Dict[str, Any]]:
    df = df_subjects.copy()
    df["accuracy"] = (df["correct"] / df["attempted"]).round(3) * 100
    df["avg_time_per_q"] = (df["time_taken"] / df["attempted"]).round(1)
    return df[
        [
            "subject",
            "accuracy",
            "marks_scored",
            "attempted",
            "correct",
            "avg_time_per_q",
            "time_taken",
        ]
    ].to_dict(orient="records")


def summarize_chapters(df_questions: pd.DataFrame) -> List[Dict[str, Any]]:
    g = df_questions.groupby(["subject", "chapter"])
    summary = g.agg(
        attempted=("correct", "size"),
        correct=("correct", "sum"),
        time_taken=("time_taken", "sum"),
    ).reset_index()
    summary["accuracy"] = (summary["correct"] / summary["attempted"]).round(3) * 100
    summary["avg_time_per_q"] = (summary["time_taken"] / summary["attempted"]).round(1)

    # difficulty distribution
    diff_counts = (
        df_questions.groupby(["subject", "chapter", "difficulty"])
        .size()
        .unstack(fill_value=0)
    )
    summary = summary.join(diff_counts, on=["subject", "chapter"]).fillna(0)
    return summary.to_dict(orient="records")


def prepare_llm_context(raw: Dict[str, Any]) -> Dict[str, Any]:
    """High-level orchestrator – returns dict ready for GPT JSON-mode."""
    df_subj, df_q = build_dataframes(raw)

    context = {
        "overall": {
            "score": raw["totalMarkScored"],
            "max_score": raw["test"]["totalMarks"],
            "accuracy": round(raw["accuracy"], 2),
            "time_taken_sec": raw["totalTimeTaken"],
            "total_questions": raw["test"]["totalQuestions"],
        },
        "subjects": summarize_subjects(df_subj),
        "chapters": summarize_chapters(df_q),
        # quick pointers for the prompt…
        "strongest_subject": max(df_subj.to_dict("records"), key=lambda d: d["accuracy"])[
            "subject"
        ],
        "weakest_chapter": min(
            summarize_chapters(df_q), key=lambda c: c["accuracy"]
        )["chapter"],
    }
    return context


# ---------------------------------------------------------------------------
# CLI entry-point (handy during dev)
# ---------------------------------------------------------------------------

def main(path: str | Path) -> None:
    raw = json.loads(Path(path).read_text())
    ctx = prepare_llm_context(raw)
    print(json.dumps(ctx, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python -m analytics <report.json>")
    main(sys.argv[1])
