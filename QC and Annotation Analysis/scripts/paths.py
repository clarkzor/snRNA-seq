from pathlib import Path

def find_project_root(start=None):
    start = Path.cwd() if start is None else Path(start)

    for path in [start] + list(start.parents):
        if (path / ".git").exists() or (path / "README.md").exists():
            return path

    raise FileNotFoundError("Could not find project root.")

PROJECT_DIR = find_project_root()
DATA_DIR = PROJECT_DIR / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
RESULTS_DIR = PROJECT_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
QC_DIR = RESULTS_DIR / "qc"
PROCESSED_DATA_DIR = RESULTS_DIR / "processed_data"
