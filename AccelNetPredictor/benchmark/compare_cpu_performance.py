#!/usr/bin/env python3
"""Compare warmed wall times and ALL energy/force/virial components.

Use executables built with the same driver/compiler/options and different
library revisions. CPU affinity, alternating order and paired ratios reduce noise;
run on an otherwise idle host. Threshold failures have a nonzero exit status.
"""
import argparse
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys


def measure(executable, family, natoms, args, mode):
    command = [str(executable), family, str(natoms), str(args.seconds), str(args.data), mode]
    result = subprocess.run(command, check=True, text=True, capture_output=True, timeout=300)
    timing, values = None, []
    for line in result.stdout.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "TIMING":
            timing = float(fields[1])
        elif fields[0] in ("ENERGY", "FORCE", "VIRIAL"):
            values.extend(float(value) for value in fields[1:])
    if timing is None or not math.isfinite(timing) or timing <= 0:
        raise ValueError("missing/invalid wall time")
    if len(values) != 1 + 3 * natoms + 9 or not all(map(math.isfinite, values)):
        raise ValueError("missing/nonfinite energy, force or virial components")
    return timing, values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().parents[1] / "test/data")
    parser.add_argument("--families", nargs="+", default=["chebyshev", "lj", "n2p2-g5", "n2p2-g4"])
    parser.add_argument("--sizes", nargs="+", type=int, default=[8, 64, 512])
    parser.add_argument("--samples", type=int, default=7)
    parser.add_argument("--seconds", type=float, default=0.2)
    parser.add_argument("--max-slowdown", type=float, default=1.10)
    parser.add_argument("--baseline-mode", choices=["structure", "reference"], default="structure")
    parser.add_argument("--candidate-mode", choices=["structure", "batch"], default="structure")
    parser.add_argument("--cpu", type=int, help="Linux CPU ID (default: first available CPU)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.samples < 3 or args.seconds <= 0 or args.max_slowdown < 1 or any(n < 1 for n in args.sizes):
        parser.error("need >=3 samples, positive duration/sizes and max-slowdown >=1")
    args.baseline = args.baseline.resolve()
    args.candidate = args.candidate.resolve()
    args.data = args.data.resolve()
    affinity = None
    if hasattr(os, "sched_getaffinity"):
        available = os.sched_getaffinity(0)
        affinity = args.cpu if args.cpu is not None else min(available)
        if affinity not in available:
            parser.error("requested CPU is outside the available affinity mask")
        os.sched_setaffinity(0, {affinity})
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    report = dict(host=platform.node(), platform=platform.platform(), cpu=affinity,
                  baseline=str(args.baseline), candidate=str(args.candidate),
                  baseline_mode=args.baseline_mode, candidate_mode=args.candidate_mode, samples=args.samples,
                  seconds_per_sample=args.seconds, max_slowdown=args.max_slowdown,
                  ratio_method="median of paired candidate/baseline wall times",
                  absolute_tolerance=2e-10, relative_tolerance=2e-10, cases=[])
    failed = False
    for family in args.families:
        for natoms in args.sizes:
            times = {"baseline": [], "candidate": []}
            reference = None
            max_error = 0.0
            correct = True
            for sample in range(args.samples):
                names = ["baseline", "candidate"] if sample % 2 == 0 else ["candidate", "baseline"]
                for name in names:
                    mode = args.candidate_mode if name == "candidate" else args.baseline_mode
                    timing, values = measure(getattr(args, name), family, natoms, args, mode)
                    times[name].append(timing)
                    if reference is None:
                        reference = values
                    for actual, expected in zip(values, reference):
                        error = abs(actual - expected)
                        max_error = max(max_error, error)
                        correct &= error <= 2e-10 + 2e-10 * abs(expected)
            # Pair adjacent measurements before taking the median. A frequency
            # transition halfway through a run can put separate medians in
            # different timing regimes, even when both binaries run equally fast.
            paired_ratios = [a / b for a, b in zip(times["candidate"], times["baseline"])]
            ratio = statistics.median(paired_ratios)
            passed = correct and ratio <= args.max_slowdown
            failed |= not passed
            case = dict(family=family, natoms=natoms, seconds=times, ratio=ratio, paired_ratios=paired_ratios,
                        max_absolute_error=max_error, correct=correct, passed=passed)
            report["cases"].append(case)
            print(f"{family:12s} N={natoms:6d} ratio={ratio:.3f} error={max_error:.3g} "
                  f"{'PASS' if passed else 'FAIL'}", flush=True)
    report["passed"] = not failed
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
