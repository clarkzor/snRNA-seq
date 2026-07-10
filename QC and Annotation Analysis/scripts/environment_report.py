# scripts/environment_report.py

from pathlib import Path
from datetime import datetime
import sys
import platform
from importlib.metadata import version, PackageNotFoundError
import yaml


DEFAULT_PACKAGES = [
    "scanpy",
    "anndata",
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "seaborn",
    "scrublet",
    "leidenalg",
    "igraph",
    "umap-learn",
    "scikit-learn",
    "h5py",
    "pyyaml",
    "stream2",
]


def get_software_versions(packages=None):
    """
    Collect Python, OS, and package version information.
    """
    if packages is None:
        packages = DEFAULT_PACKAGES

    versions = {
        "timestamp": datetime.now().isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {},
    }

    for package in packages:
        try:
            versions["packages"][package] = version(package)
        except PackageNotFoundError:
            versions["packages"][package] = "not installed"

    return versions


def write_software_versions_txt(output_path, packages=None):
    """
    Write software versions to a human-readable text file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    versions = get_software_versions(packages)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("Software Versions\n")
        f.write("=================\n\n")
        f.write(f"Timestamp: {versions['timestamp']}\n")
        f.write(f"Python: {versions['python']}\n")
        f.write(f"Platform: {versions['platform']}\n\n")

        f.write("Packages:\n")
        for package, package_version in versions["packages"].items():
            f.write(f"- {package}: {package_version}\n")

    print(f"Software versions written to: {output_path}")


def write_run_config_yaml(output_path, run_info=None, parameters=None, packages=None):
    """
    Write run-specific configuration information to a YAML file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    versions = get_software_versions(packages)

    run_config = {
        "run_metadata": {
            "timestamp": versions["timestamp"],
            "python": versions["python"],
            "platform": versions["platform"],
        },
        "software": versions["packages"],
        "run_info": run_info or {},
        "parameters": parameters or {},
    }

    with open(output_path, "w", encoding="utf-8") as f:
        yaml.dump(run_config, f, sort_keys=False)

    print(f"Run configuration written to: {output_path}")
