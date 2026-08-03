#!/usr/bin/env python3
"""Summarize Linux LAMMPS benchmark logs as median and IQR tables."""

from __future__ import annotations

import argparse
import re
import statistics
from pathlib import Path


NAME_PATTERN = re.compile(
    r"(?P<system>water|tio2)-(?P<size>small|medium|large)-(?P<atoms>\d+)-"
    r"(?P<mode>fixed|dynamic)-(?P<pair>aenet|accelnet)-r(?P<ranks>\d+)-t(?P<trial>\d+)\.log"
)


def percentile(values, fraction):
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def parse_log(path: Path):
    match = NAME_PATTERN.fullmatch(path.name)
    if not match:
        return None
    text = path.read_text(encoding="utf-8")
    loops = re.findall(
        r"Loop time of\s+([0-9.eE+-]+)\s+on\s+(\d+)\s+procs for\s+(\d+)\s+steps with\s+(\d+)\s+atoms",
        text,
    )
    if not loops:
        raise ValueError(f"measured loop time not found in {path}")
    loop, ranks, steps, atoms = loops[-1]
    metadata = match.groupdict()
    implementation = metadata.pop("pair")
    record = {
        **metadata,
        "implementation": implementation,
        "path": str(path),
        "loop": float(loop),
        "steps": int(steps),
    }
    if int(ranks) != int(record["ranks"]) or int(atoms) != int(record["atoms"]):
        raise ValueError(f"filename metadata differs from LAMMPS output: {path}")
    for section in ("Pair", "Neigh", "Comm", "Modify", "Kspace"):
        rows = re.findall(
            rf"^{section}\s*\|\s*[^|]+\|\s*([0-9.eE+-]+)\s*\|",
            text,
            re.MULTILINE,
        )
        record[f"{section.lower()}_time"] = float(rows[-1]) if rows else 0.0
    other_rows = re.findall(
        r"^Other\s*\|\s*\|\s*([0-9.eE+-]+)\s*\|",
        text,
        re.MULTILINE,
    )
    record["other_time"] = float(other_rows[-1]) if other_rows else 0.0
    dangerous_rows = re.findall(r"^Dangerous builds\s*=\s*(\d+)", text, re.MULTILINE)
    record["dangerous_builds"] = int(dangerous_rows[-1]) if dangerous_rows else -1
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("log_directory", type=Path)
    args = parser.parse_args()
    records = [parse_log(path) for path in sorted(args.log_directory.glob("*.log"))]
    records = [record for record in records if record is not None]
    if not records:
        raise SystemExit("no benchmark logs found")

    groups = {}
    for record in records:
        key = (
            record["system"], record["size"], int(record["atoms"]), record["mode"],
            record["implementation"], int(record["ranks"]), record["steps"],
        )
        groups.setdefault(key, []).append(record)

    print("system\tsize\tatoms\tmode\tpair\tranks\tsteps\ttrials\tloop_median_s\tloop_q1_s\tloop_q3_s\tsteps_per_s\tatom_steps_per_s\tpair_median_s\tneigh_median_s\tcomm_median_s\tmodify_median_s\tkspace_median_s\tother_median_s\tmax_dangerous_builds")
    for key in sorted(groups, key=lambda item: (item[0], item[2], item[3], item[5], item[4])):
        group = groups[key]
        loops = [record["loop"] for record in group]
        median_loop = statistics.median(loops)
        steps_per_second = key[6] / median_loop
        timing = [statistics.median(record[f"{section}_time"] for record in group)
                  for section in ("pair", "neigh", "comm", "modify", "kspace", "other")]
        print(
            f"{key[0]}\t{key[1]}\t{key[2]}\t{key[3]}\t{key[4]}\t{key[5]}\t{key[6]}\t{len(group)}\t"
            f"{median_loop:.9g}\t{percentile(loops, 0.25):.9g}\t{percentile(loops, 0.75):.9g}\t"
            f"{steps_per_second:.9g}\t{steps_per_second * key[2]:.9g}\t"
            + "\t".join(f"{value:.9g}" for value in timing)
            + f"\t{max(record['dangerous_builds'] for record in group)}"
        )


if __name__ == "__main__":
    main()
