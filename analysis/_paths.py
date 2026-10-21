"""Shared paths for the analysis scripts.

    analysis/scripts/*.py             → outputs in analysis/{results,figures}/
    analysis/historical/scripts/*.py  → outputs in analysis/historical/{results,figures}/

Figures (.png / .pdf) go to figures/, everything else (tables, reports,
caches) to results/. Each script imports this via:

    sys.path.insert(0, str(<analysis dir>))
    from _paths import ROOT, output
"""
from __future__ import annotations

from pathlib import Path

ANALYSIS = Path(__file__).resolve().parent
ROOT = ANALYSIS.parent
DATA = ROOT / "data"
_FIGURE_SUFFIXES = {".png", ".pdf"}


def output(name: str, historical: bool = False) -> Path:
    """Path for an analysis output (or another script's output to read)."""
    base = ANALYSIS / "historical" if historical else ANALYSIS
    sub = "figures" if Path(name).suffix in _FIGURE_SUFFIXES else "results"
    path = base / sub / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
