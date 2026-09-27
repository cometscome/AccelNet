#!/usr/bin/env python3
"""Guard G5 moment performance on the original degree-4 scaling fixture.

All four paths use one binary/model and exclude neighbor construction. Run only
with ACCELNET_TARGET_SERIAL=ON and OpenMP disabled, on an otherwise idle host.
This fixture protects the established moment advantage, not every possible model.
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--seconds", type=float, default=0.05)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--max-slowdown", type=float, default=1.10)
    args = parser.parse_args()
    if args.seconds <= 0 or args.rounds < 2 or args.max_slowdown < 1:
        parser.error("need positive seconds, at least two rounds and max-slowdown >= 1")
    exe = args.executable.resolve()
    symbols = subprocess.check_output(["nm", "-u", str(exe)], text=True)
    if re.search(r"\b(?:GOMP_|__kmpc_|__nvomp_|omp_)[A-Za-z0-9_]*", symbols):
        parser.error("use a binary built without OpenMP, not just OMP_NUM_THREADS=1")
    cpu = None
    if hasattr(os, "sched_getaffinity"):
        available = os.sched_getaffinity(0)
        cpu = args.cpu if args.cpu is not None else min(available)
        if cpu not in available:
            parser.error("requested CPU is outside the available affinity mask")
        os.sched_setaffinity(0, {cpu})
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    rows = []
    for round_id in range(args.rounds):
        for mode in ([1, 3] if round_id % 2 == 0 else [3, 1]):
            command = [str(exe), "64", "4", str(args.seconds), str(mode), "1.7",
                       "cpu-shared", "g5-scaling", "no-neighbors", "64"]
            result = subprocess.run(command, env=env, check=True, capture_output=True,
                                    text=True, timeout=300)
            samples = [list(map(float, line.split()[2:4]))
                       for line in result.stdout.splitlines() if line.startswith("TIMING ")]
            errors = [float(line.split()[4]) for line in result.stdout.splitlines()
                      if line.startswith("CASE ")]
            if len(samples) != 5 or len(errors) != 1 or not all(
                    math.isfinite(t) and t > 0 for pair in samples for t in pair):
                raise ValueError("missing/invalid benchmark results")
            if not math.isfinite(errors[0]) or errors[0] > 2e-9:
                raise ValueError("energy/force/virial mismatch")
            rows.append(dict(round=round_id, mode=mode, command=command, samples=samples,
                             error=errors[0], old=statistics.median(s[0] for s in samples),
                             common=statistics.median(s[1] for s in samples)))
    direct = [r for r in rows if r["mode"] == 1]
    moment = [r for r in rows if r["mode"] == 3]
    ratios = dict(
        old_moment_over_direct=statistics.median(m["old"] / d["old"] for d, m in zip(direct, moment)),
        common_moment_over_old_moment=statistics.median(m["common"] / m["old"] for m in moment),
        common_moment_over_common_direct=statistics.median(
            m["common"] / d["common"] for d, m in zip(direct, moment)))
    passed = (ratios["common_moment_over_old_moment"] <= args.max_slowdown
              and ratios["common_moment_over_common_direct"] <= 1.0)
    report = dict(cpu=cpu, rows=rows, ratios=ratios, max_slowdown=args.max_slowdown, passed=passed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(ratios, indent=2))
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
