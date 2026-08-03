#!/usr/bin/env python3
"""Compare run-0 energy and force dumps from two LAMMPS pair styles."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


def read_energy(path: Path) -> float:
    lines = path.read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(lines):
        if re.fullmatch(r"\s*Step\s+PotEng\s*", line):
            return float(lines[index + 1].split()[1])
    raise ValueError(f"thermodynamic energy not found in {path}")


def read_forces(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    forces = {}
    for line in lines[header + 1:]:
        fields = line.split()
        if len(fields) == 8:
            forces[int(fields[0])] = tuple(float(value) for value in fields[5:8])
    return forces


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_log", type=Path)
    parser.add_argument("reference_dump", type=Path)
    parser.add_argument("candidate_log", type=Path)
    parser.add_argument("candidate_dump", type=Path)
    parser.add_argument("--tolerance", type=float, default=1.0e-8)
    args = parser.parse_args()

    reference_energy = read_energy(args.reference_log)
    candidate_energy = read_energy(args.candidate_log)
    reference_forces = read_forces(args.reference_dump)
    candidate_forces = read_forces(args.candidate_dump)
    if reference_forces.keys() != candidate_forces.keys():
        raise SystemExit("FAIL: atom IDs differ")
    energy_error = abs(reference_energy - candidate_energy)
    force_error = max(
        abs(left - right)
        for atom_id in reference_forces
        for left, right in zip(reference_forces[atom_id], candidate_forces[atom_id])
    )
    print(f"REFERENCE_ENERGY {reference_energy:.15g}")
    print(f"CANDIDATE_ENERGY {candidate_energy:.15g}")
    print(f"MAX_ABS_ENERGY_ERROR {energy_error:.6e}")
    print(f"MAX_ABS_FORCE_ERROR {force_error:.6e}")
    if not math.isfinite(energy_error + force_error):
        raise SystemExit("FAIL: non-finite comparison result")
    if energy_error > args.tolerance or force_error > args.tolerance:
        raise SystemExit(f"FAIL: difference exceeds tolerance {args.tolerance:g}")
    print("PASS")


if __name__ == "__main__":
    main()

