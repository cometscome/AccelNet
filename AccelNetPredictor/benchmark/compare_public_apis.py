#!/usr/bin/env python3
"""Compare the identical public-API driver linked to old/new serial libraries."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ('before', 'after', 'aenet', 'n2p2', 'output'):
        p.add_argument('--'+k, type=Path, required=True)
    p.add_argument('--core', type=int, default=6)
    p.add_argument('--rounds', type=int, default=4)
    p.add_argument('--seconds', type=float, default=.15)
    p.add_argument('--aenet-energy-only', action='store_true',
                   help='check both energy-only APIs at multiple atom counts')
    p.add_argument('--max-energy-ratio', type=float, default=None,
                   help='optional after/before limit for aenet energy-only calls')
    a = p.parse_args()
    if a.max_energy_ratio is not None and a.max_energy_ratio <= 0: p.error('ratio must be positive')
    if a.rounds < 2 or a.seconds <= 0: p.error('need >=2 rounds and positive seconds')
    a.output.mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, {a.core})
    binaries = {'before': a.before.resolve(), 'after': a.after.resolve()}
    for exe in binaries.values():
        undefined = subprocess.check_output(['nm', '-u', exe], text=True)
        if any(x in undefined for x in ('GOMP_', '__kmpc_', '__nvomp')):
            raise RuntimeError('CPU baseline must compile OpenMP out')
    report = dict(core=a.core, rounds=a.rounds, seconds=a.seconds,
                  binaries={k: dict(path=str(v), sha256=hashlib.sha256(v.read_bytes()).hexdigest())
                            for k, v in binaries.items()}, rows=[])
    reference = {}
    cases = [('aenet', a.aenet, 192), ('n2p2', a.n2p2, 512),
             ('combined', Path('.'), 128), ('multi-chebyshev', Path('.'), 128),
             ('mixed-components', Path('.'), 128)]
    if a.aenet_energy_only: cases = [('aenet', a.aenet, n) for n in (64, 192, 512)]
    for family, directory, n in cases:
        paths = ['structure', 'energy', 'batch']
        if family in ('aenet', 'n2p2'): paths += ['atomic', 'atomic-energy']
        if a.aenet_energy_only: paths = ['energy', 'atomic-energy']
        for api in paths:
            row = dict(family=family, atoms=n, api=api, samples={k: [] for k in binaries}, runs=[])
            for repeat in range(a.rounds):
                for label in (['before', 'after'] if repeat % 2 == 0 else ['after', 'before']):
                    cmd = [str(binaries[label]), family, str(directory.resolve()), str(n), str(a.seconds), api]
                    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
                    result = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)
                    (a.output/f'{family}-{n}-{api}-{repeat}-{label}.log').write_text(result.stdout+result.stderr)
                    if result.returncode: raise RuntimeError(result.stdout[-2000:]+result.stderr[-2000:])
                    lines = result.stdout.splitlines()
                    timing = [float(l.split()[1]) for l in lines if l.startswith('TIMING ')]
                    energy = [float(l.split()[1]) for l in lines if l.startswith('ENERGY ')]
                    fw = [float(x) for l in lines if l.startswith(('FORCE ', 'VIRIAL ')) for x in l.split()[1:]]
                    if len(timing) != 1 or len(energy) != 1 or len(fw) != 3*n+9: raise RuntimeError('incomplete output')
                    values = energy if 'energy' in api else energy+fw
                    key = (family, n, 'energy' if 'energy' in api else 'force')
                    expected = reference.setdefault(key, values)
                    errors = [abs(x-y)/(2e-10+2e-10*abs(y)) for x,y in zip(values, expected)]
                    if not all(e <= 1 for e in errors): raise AssertionError((family, api, label, max(errors)))
                    # Also compare every energy-only result to the force path's energy.
                    if (family, n, 'force') in reference and abs(energy[0]-reference[(family,n,'force')][0]) > 2e-10*(1+abs(energy[0])):
                        raise AssertionError('energy-only differs')
                    row['samples'][label].append(timing[0])
                    row['runs'].append(dict(command=cmd, label=label, max_tolerance_fraction=max(errors)))
            row['seconds'] = {k: statistics.median(v) for k,v in row['samples'].items()}
            row['after_over_before'] = row['seconds']['after']/row['seconds']['before']
            report['rows'].append(row)
            if family == 'aenet' and 'energy' in api and a.max_energy_ratio is not None:
                row['performance_passed'] = row['after_over_before'] <= a.max_energy_ratio
            print(family, n, api, row['seconds'], row['after_over_before'], flush=True)
            (a.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    report['numerical_passed'] = True
    report['max_energy_ratio'] = a.max_energy_ratio
    report['passed'] = all(r.get('performance_passed', True) for r in report['rows'])
    (a.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    if not report['passed']: raise SystemExit('energy performance regression; see report.json')


if __name__ == '__main__':
    main()
