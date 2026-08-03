# Repository validation outputs

This directory stores validation reports for the repository's Jupyter notebook
files. Regenerate them with:

```bash
python scripts/validate_notebooks.py
```

The validator checks:

1. Jupyter notebook schema and format using `nbformat.validate`.
2. Static Python syntax for every code cell after IPython transformations.

It does not execute notebook cells or validate data-dependent Scanpy results.
The stage-specific `results/<stage>[/<condition>]/validation/` directories are
reserved for biological and runtime validation outputs.
