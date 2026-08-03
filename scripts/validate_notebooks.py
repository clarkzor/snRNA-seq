"""Validate repository notebooks and write reproducible reports.

This script performs two static checks for every ``notebooks/*.ipynb`` file:

1. Jupyter notebook-format/schema validation with ``nbformat.validate``.
2. Python syntax validation for each code cell after IPython input transforms.

It does not execute the notebooks and therefore does not validate Scanpy runtime
behavior, data availability, or biological results.
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import nbformat
from IPython.core.inputtransformer2 import TransformerManager


def find_project_root(start: Path) -> Path:
    """Find the repository root using project-specific marker files."""
    markers = (
        Path("README.md"),
        Path("notebooks"),
        Path("scripts") / "output_utils.py",
    )

    for candidate in (start.resolve(), *start.resolve().parents):
        if all((candidate / marker).exists() for marker in markers):
            return candidate

    raise FileNotFoundError(
        "Could not locate the repository root. Expected README.md, notebooks/, "
        "and scripts/output_utils.py in the same directory."
    )


def validate_notebook(path: Path, transformer: TransformerManager) -> dict[str, Any]:
    """Validate one notebook and return a serializable result record."""
    result: dict[str, Any] = {
        "notebook": path.name,
        "relative_path": str(path),
        "nbformat_valid": False,
        "nbformat_error": "",
        "code_cells": 0,
        "python_syntax_valid": False,
        "syntax_error_count": 0,
        "syntax_errors": [],
        "status": "FAIL",
    }

    try:
        notebook = nbformat.read(path, as_version=4)
        nbformat.validate(notebook)
        result["nbformat_valid"] = True
    except Exception as exc:  # schema/read failures need to be captured in report
        result["nbformat_error"] = f"{type(exc).__name__}: {exc}"
        return result

    syntax_errors: list[dict[str, Any]] = []

    for cell_index, cell in enumerate(notebook.cells):
        if cell.cell_type != "code":
            continue

        result["code_cells"] += 1
        source = cell.source or ""

        try:
            transformed = transformer.transform_cell(source)
            compile(transformed, f"{path.name}:cell_{cell_index}", "exec")
        except SyntaxError as exc:
            syntax_errors.append(
                {
                    "cell_index": cell_index,
                    "line": exc.lineno,
                    "offset": exc.offset,
                    "message": exc.msg,
                    "text": (exc.text or "").strip(),
                }
            )
        except Exception as exc:
            syntax_errors.append(
                {
                    "cell_index": cell_index,
                    "line": None,
                    "offset": None,
                    "message": f"{type(exc).__name__}: {exc}",
                    "text": "",
                }
            )

    result["syntax_errors"] = syntax_errors
    result["syntax_error_count"] = len(syntax_errors)
    result["python_syntax_valid"] = not syntax_errors
    result["status"] = (
        "PASS"
        if result["nbformat_valid"] and result["python_syntax_valid"]
        else "FAIL"
    )
    return result


def write_reports(
    results: list[dict[str, Any]],
    output_dir: Path,
    project_root: Path,
) -> None:
    """Write CSV, JSON, and readable text validation reports."""
    output_dir.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now(timezone.utc).isoformat()

    csv_path = output_dir / "notebook_format_validation.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "notebook",
                "relative_path",
                "nbformat_valid",
                "code_cells",
                "python_syntax_valid",
                "syntax_error_count",
                "status",
                "nbformat_error",
            ],
        )
        writer.writeheader()
        for result in results:
            row = {key: result.get(key, "") for key in writer.fieldnames}
            writer.writerow(row)

    json_path = output_dir / "notebook_format_validation.json"
    json_payload = {
        "generated_at_utc": generated_at,
        "project_root": ".",
        "checks": {
            "notebook_schema": "nbformat.validate",
            "code_cell_syntax": "IPython transform followed by Python compile",
            "notebook_execution": False,
        },
        "summary": {
            "notebooks_checked": len(results),
            "passed": sum(result["status"] == "PASS" for result in results),
            "failed": sum(result["status"] != "PASS" for result in results),
        },
        "results": results,
    }
    json_path.write_text(json.dumps(json_payload, indent=2), encoding="utf-8")

    text_path = output_dir / "notebook_format_validation.txt"
    lines = [
        "Notebook format and static syntax validation",
        "=" * 45,
        f"Generated (UTC): {generated_at}",
        "Repository: .",
        "",
        "Checks performed:",
        "- Jupyter notebook schema validation with nbformat.validate",
        "- Static Python syntax compilation of every code cell after IPython transforms",
        "- Notebook execution: NOT performed",
        "",
    ]

    for result in results:
        lines.append(
            f"[{result['status']}] {result['notebook']} | "
            f"nbformat={result['nbformat_valid']} | "
            f"code_cells={result['code_cells']} | "
            f"syntax={result['python_syntax_valid']}"
        )
        if result["nbformat_error"]:
            lines.append(f"  nbformat error: {result['nbformat_error']}")
        for error in result["syntax_errors"]:
            lines.append(
                "  syntax error: "
                f"cell={error['cell_index']}, line={error['line']}, "
                f"message={error['message']}"
            )

    passed = sum(result["status"] == "PASS" for result in results)
    failed = len(results) - passed
    lines.extend(
        [
            "",
            f"Summary: {passed} passed; {failed} failed; {len(results)} checked.",
            "",
            "Important: PASS confirms valid notebook structure and static code syntax only. ",
            "It does not confirm that cells execute successfully against the required .h5ad files.",
        ]
    )
    text_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate all Jupyter notebooks in this repository."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=None,
        help="Repository root. Defaults to automatic detection.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Report directory. Defaults to <project-root>/validation.",
    )
    args = parser.parse_args()

    project_root = (
        args.project_root.resolve()
        if args.project_root is not None
        else find_project_root(Path.cwd())
    )
    notebooks_dir = project_root / "notebooks"
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else project_root / "validation"
    )

    notebook_paths = sorted(notebooks_dir.glob("*.ipynb"))
    if not notebook_paths:
        raise FileNotFoundError(f"No .ipynb files found in {notebooks_dir}")

    transformer = TransformerManager()
    results = []
    for path in notebook_paths:
        result = validate_notebook(path, transformer)
        result["relative_path"] = str(path.relative_to(project_root))
        results.append(result)

    write_reports(results, output_dir, project_root)

    for result in results:
        print(f"{result['status']:4}  {result['relative_path']}")

    failures = [result for result in results if result["status"] != "PASS"]
    print(f"\nReports written to: {output_dir}")
    print(f"Summary: {len(results) - len(failures)} passed; {len(failures)} failed.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
