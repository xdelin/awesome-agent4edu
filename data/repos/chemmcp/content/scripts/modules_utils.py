import inspect
import subprocess
from datetime import datetime, timezone, timedelta

def get_class_file_path(cls) -> str:
    """
    Return the source file path where the given class is defined.
    Raises ValueError if the file cannot be determined.
    """
    path = inspect.getsourcefile(cls) or inspect.getfile(cls)
    if path is None:
        raise ValueError(f"Cannot locate source file for class {cls.__name__}")
    return path

def git_last_commit_date_for_class(cls) -> str:
    """
    Return the date of the last Git commit for the file defining cls,
    formatted as 'yyyy/MM/DD' in UTC−5. If the date cannot be determined,
    return 'Unknown'.
    """
    try:
        file_path = get_class_file_path(cls)
    except ValueError:
        return "Unknown"

    try:
        # Use %cI for ISO 8601 commit date (e.g. 2025-05-10T18:23:45+02:00)
        iso_ts = subprocess.check_output(
            ['git', 'log', '-1', '--format=%cI', '--', file_path],
            stderr=subprocess.DEVNULL,
            text=True
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        # CalledProcessError: git returns non-zero; FileNotFoundError: git not installed
        return "Unknown"

    try:
        # Parse the ISO timestamp into a datetime with tzinfo
        dt = datetime.fromisoformat(iso_ts)
        # Convert to UTC−5
        target_tz = timezone(timedelta(hours=-5))
        dt_utc5 = dt.astimezone(target_tz)
        # Format as yyyy/MM/DD
        return dt_utc5.strftime('%Y/%m/%d')
    except Exception:
        return "Unknown"

