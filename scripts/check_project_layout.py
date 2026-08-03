"""Validate the repository layout before running analysis notebooks."""

from __future__ import annotations

from pathlib import Path
import sys


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]

    required = [
        project_root / "README.md",
        project_root / "notebooks",
        project_root / "scripts" / "output_utils.py",
        project_root / "data",
        project_root / "results",
    ]

    missing = [path for path in required if not path.exists()]
    if missing:
        print("Repository layout check failed. Missing:")
        for path in missing:
            print(f"  - {path}")
        return 1

    legacy_nested = [
        project_root / "notebooks" / "data",
        project_root / "notebooks" / "results",
    ]
    present_legacy = [path for path in legacy_nested if path.exists()]

    print("Repository root:", project_root)
    print("Top-level data directory:", project_root / "data")
    print("Top-level results directory:", project_root / "results")

    if present_legacy:
        print("\nERROR: legacy nested directories were detected:")
        for path in present_legacy:
            print(f"  - {path}")
        print(
            "\nMove their contents to the corresponding top-level data/ or "
            "results/ directory, then remove the nested directories."
        )
        return 2

    print("\nPASS: no notebooks/data or notebooks/results directories exist.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
