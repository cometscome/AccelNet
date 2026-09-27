#!/usr/bin/env python3
"""Compare saved descriptor binaries in alternating order, with full E/F/W checks.

Example: --variant before host /path/to/before --variant after host /path/to/after
Use serial builds for CPU comparisons. GPU selection is inherited from the
environment. Each binary performs five alternating reference/candidate samples.
Chebyshev is the default; use --family lj/g4/g5/behler --modes 1 for generic
direct kernels, or --modes 0 to include the CPU reference's automatic G5 policy.
"""
import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variant', nargs=3, action='append', required=True,
                        metavar=('LABEL', 'BACKEND', 'EXECUTABLE'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sizes', type=int, nargs='+', default=[512, 4096])
    parser.add_argument('--orders', type=int, nargs='+', default=[5])
    parser.add_argument('--modes', type=int, nargs='+', default=[1, 2])
    parser.add_argument('--family', choices=['chebyshev', 'lj', 'g4', 'g4-distinct', 'g5', 'behler', 'lj-behler', 'g4-series', 'g5-series', 'g5-scaling'],
                        default='chebyshev')
    parser.add_argument('--neighbors', type=int, default=0, help='fixed original G5 scaling environment size')
    parser.add_argument('--rounds', type=int, default=3)
    parser.add_argument('--seconds', type=float, default=0.12)
    parser.add_argument('--cpu', type=int, default=6)
    args = parser.parse_args()
    if args.rounds < 1 or args.seconds <= 0:
        parser.error('rounds and seconds must be positive')
    if args.family != 'chebyshev' and any(mode not in (0, 1, 2, 3) for mode in args.modes):
        parser.error('generic modes: 0 auto, 1 direct, 2 moments with threshold, 3 forced moments')
    labels = [v[0] for v in args.variant]
    if len(set(labels)) != len(labels):
        parser.error('variant labels must be unique')
    args.output.mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, {args.cpu})
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               OMP_TARGET_OFFLOAD='MANDATORY')
    report = []
    for n in args.sizes:
        for order in args.orders:
            for mode in args.modes:
                for repeat in range(args.rounds):
                    variants = args.variant if repeat % 2 == 0 else args.variant[::-1]
                    for label, backend, exe in variants:
                        command = [str(Path(exe).resolve()), str(n), str(order), str(args.seconds),
                                   str(mode), '1.7', backend, args.family, 'no-neighbors']
                        if args.neighbors:
                            command.append(str(args.neighbors))
                        result = subprocess.run(command, env=env, text=True, capture_output=True,
                                                timeout=600)
                        name = f'{label}-n{n}-order{order}-mode{mode}-round{repeat}'
                        if args.family != 'chebyshev':
                            name = args.family + '-' + name
                        (args.output / f'{name}.log').write_text(result.stdout + result.stderr)
                        if result.returncode:
                            raise RuntimeError(f'{name} failed; see log')
                        samples = [list(map(float, line.split()[2:]))
                                   for line in result.stdout.splitlines() if line.startswith('TIMING ')]
                        if len(samples) != 5:
                            raise RuntimeError(f'{name}: missing five timing samples')
                        case = next(line.split() for line in result.stdout.splitlines()
                                    if line.startswith('CASE '))
                        phases = [list(map(float, line.split()[3:]))
                                  for line in result.stdout.splitlines() if line.startswith('PROFILE ')]
                        item = dict(label=label, backend=backend, family=args.family, natoms=n, order=order, mode=mode,
                                    repeat=repeat, command=command, neighbors=int(case[3]),
                                    max_error=float(case[4]), samples=samples,
                                    seconds=statistics.median(s[1] for s in samples),
                                    reference_seconds=statistics.median(s[0] for s in samples),
                                    ratio_to_reference=statistics.median(s[1]/s[0] for s in samples),
                                    phase_medians=[statistics.median(p[i] for p in phases)
                                                   for i in range(8)])
                        report.append(item)
                        (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
                        print(name, 'seconds', item['seconds'], 'ratio', item['ratio_to_reference'],
                              flush=True)


if __name__ == '__main__':
    main()
