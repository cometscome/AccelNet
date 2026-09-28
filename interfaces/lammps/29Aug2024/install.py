#!/usr/bin/env python3
"""Compatibility entry point for the shared LAMMPS installer."""
import runpy
from pathlib import Path

if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parent.parent / "install.py"),
                   run_name="__main__")
