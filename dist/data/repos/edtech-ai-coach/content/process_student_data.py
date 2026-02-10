from __future__ import annotations
from typing import List, Dict, Any
from collections import defaultdict, Counter
import statistics
from bs4 import BeautifulSoup

def _extract_subject_names_from_syllabus(html: str) -> List[str]:
    soup = BeautifulSoup(html, "html.parser")
    return [h2.get_text(strip=True) for h2 in soup.find_all("h2")]

def _build_subject_id_to_name_map(subjects: List[Dict[str, Any]], names: List[str]) -> Dict[str, str]:
    mapping = {}
    for idx, subj in enumerate(subjects):
        oid = subj["subjectId"]["$oid"]
        mapping[oid] = names[idx] if idx < len(names) else f"Subject_{idx+1}"
    return mapping

def _calculate_performance_tier(accuracy: float) -> str:
    if accuracy >= 80:
        return "strong"
    elif 60 <= accuracy < 80:
        return "average"
    return "needs_improvement"

def _generate_actionable_suggestions(subjects: List[Dict], chapters: List[Dict], overall_acc: float) -> List[str]:
    tier = _calculate_performance_tier(overall_acc)
    suggestions = []

    accuracy_values = [c["accuracy"] for c in chapters if c.get("attempted", 0) > 0]
    weak_threshold = statistics.quantiles(accuracy_values, n=4)[0] if accuracy_values else 50
    weak_subjs = [s["name"] for s in subjects if s["accuracy"] < weak_threshold]

    if weak_subjs:
        suggestions.append(
            f"Focus on **{', '.join(weak_subjs[:3])}** through: "
            f"{'concept reviews' if tier == 'needs_improvement' else 'advanced problems'}"
        )
    else:
        suggestions.append("Maintain strength across all subjects with varied practice sets")

    time_suggestions = {
        "strong": "Challenge yourself with timed practice tests",
        "average": "Use a timer during practice to optimize pacing",
        "needs_improvement": "Focus on accuracy first, then gradually reduce time per question"
    }
    strategy_map = {
        "strong": "Analyze complex multi-step problems to push your limits",
        "average": "Mix practice between strong and weak areas daily",
        "needs_improvement": "Start with foundational concepts and progress gradually"
    }

    suggestions.append(time_suggestions[tier])
    suggestions.append(strategy_map[tier])
    return suggestions[:3]

def process_reports_to_summary(reports: List[Dict[str, Any]]) -> Dict[str, Any]:
    total_score = 0
    max_score = 0
    subj_stats = defaultdict(lambda: Counter(score=0, attempted=0, correct=0, time_taken=0))
    chap_stats = defaultdict(lambda: Counter(score=0, attempted=0, correct=0, time_taken=0))
    chap_levels = defaultdict(Counter)
    q_times: List[int] = []
    q_correct: List[int] = []

    for rpt in reports:
        total_score += rpt.get("totalMarkScored", 0)
        max_score += rpt.get("test", {}).get("totalMarks", 0)

        subj_names = _extract_subject_names_from_syllabus(rpt["test"].get("syllabus", ""))
        id2name = _build_subject_id_to_name_map(rpt.get("subjects", []), subj_names)

        for s in rpt.get("subjects", []):
            name = id2name.get(s["subjectId"]["$oid"], "Unknown")
            subj_stats[name]["score"] += s.get("totalMarkScored", 0)
            subj_stats[name]["attempted"] += s.get("totalAttempted", 0)
            subj_stats[name]["correct"] += s.get("totalCorrect", 0)
            subj_stats[name]["time_taken"] += s.get("totalTimeTaken", 0)

        for section in rpt.get("sections", []):
            subject_from_section = section["sectionId"]["title"].split()[0]
            for q in section.get("questions", []):
                q_meta = q["questionId"]
                chapter = q_meta["chapters"][0]["title"] if q_meta.get("chapters") else "Unknown Chapter"
                chapter_key = f"{subject_from_section} | {chapter}"

                level = q_meta.get("level", "medium")
                chap_levels[chapter_key][level] += 1

                is_correct = (
                    isinstance(q.get("markedOptions"), list) and
                    len(q["markedOptions"]) > 0 and
                    q["markedOptions"][0].get("isCorrect", False)
                )

                chap_stats[chapter_key]["attempted"] += 1
                chap_stats[chapter_key]["correct"] += int(is_correct)
                chap_stats[chapter_key]["time_taken"] += q.get("timeTaken", 0)
                chap_stats[chapter_key]["score"] += 4 if is_correct else 0

                q_times.append(q.get("timeTaken", 0))
                q_correct.append(1 if is_correct else 0)

    subjects = []
    for name, c in subj_stats.items():
        attempted = c["attempted"] or 1
        accuracy = round(c["correct"] / attempted * 100, 2)
        subjects.append({
            "name": name,
            "accuracy": accuracy,
            "time_taken": c["time_taken"],
            "score": c["score"],
            "attempted": c["attempted"],
            "correct": c["correct"]
        })

    chapters = []
    for ck, c in chap_stats.items():
        attempted = c["attempted"] or 1
        accuracy = round(c["correct"] / attempted * 100, 2)
        diff_mode = max(chap_levels[ck], key=chap_levels[ck].get) if chap_levels[ck] else "medium"
        chapters.append({
            "name": ck,
            "accuracy": accuracy,
            "time_taken": c["time_taken"],
            "difficulty": diff_mode,
            "attempted": c["attempted"]
        })

    overall_acc = round(total_score / max_score * 100, 2) if max_score else 0.0

    avg_time = round(sum(q_times) / len(q_times), 2) if q_times else 0
    avg_correct = round(sum(q_correct) / len(q_correct) * 100, 2) if q_correct else 0

    sample_times = q_times[:20]
    sample_correct = q_correct[:20]

    return {
        "overall_accuracy": overall_acc,
        "total_score": total_score,
        "max_score": max_score,
        "subjects": subjects,
        "chapters": chapters,
        "time_accuracy_raw": {
            "avg_time_per_question": avg_time,
            "avg_accuracy_percent": avg_correct,
            "sample_times_sec": sample_times,
            "sample_correctness": sample_correct
        },
        "actionable_suggestions": _generate_actionable_suggestions(subjects, chapters, overall_acc),
        "performance_tier": _calculate_performance_tier(overall_acc)
    }
