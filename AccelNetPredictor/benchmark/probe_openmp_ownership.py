#!/usr/bin/env python3
"""Diagnose host ownership schedules with clocks and live worker counts.

Linux-only supplemental probe; the monitor runs on a separate core. This does
not change the governor or hardware clocks. All numerical checks remain in the
target benchmark executable. Include the optional persistent-pool trial to
compare plain host parallel-do against the host target-team schedule.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('before-serial', 'serial', 'before-openmp', 'openmp', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--persistent-openmp', type=Path)
    p.add_argument('--cores', nargs='+', type=int, default=list(range(6, 14)))
    p.add_argument('--monitor-core', type=int, default=15)
    p.add_argument('--rounds', type=int, default=4)
    a = p.parse_args()
    if len(set(a.cores)) < 8 or a.monitor_core in a.cores or a.rounds < 2:
        p.error('need eight distinct worker cores, a separate monitor core and >=2 rounds')
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, {a.monitor_core})
    paths = [('before_off1', a.before_serial, 1), ('off1', a.serial, 1),
             ('before_on8', a.before_openmp, 8), ('on8', a.openmp, 8)]
    if a.persistent_openmp: paths.append(('persistent_on8', a.persistent_openmp, 8))
    records = {k: dict(binary=str(exe.resolve()), sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                       runs=[], samples=[], errors=[]) for k, exe, _ in paths}
    report = dict(cores=a.cores, monitor_core=a.monitor_core, rounds=a.rounds,
                  case=dict(family='g5-series', atoms=2048, order=4, mode=1, neighbors=64), records=records)
    for repeat in range(a.rounds):
        for key, exe, nt in (paths if repeat % 2 == 0 else reversed(paths)):
            env = {k: v for k, v in os.environ.items() if not k.startswith(('OMP_', 'GOMP_'))}
            cores = a.cores[:nt]
            controls = dict(OMP_NUM_THREADS=str(nt), OMP_DYNAMIC='FALSE', OMP_PROC_BIND='close',
                            OMP_PLACES=','.join('{'+str(c)+'}' for c in cores), OMP_MAX_ACTIVE_LEVELS='1',
                            OMP_WAIT_POLICY='PASSIVE', GOMP_SPINCOUNT='300000', OMP_TARGET_OFFLOAD='DISABLED',
                            OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
            env.update(controls)
            cmd = ['taskset', '-c', ','.join(map(str, cores)), str(exe.resolve()),
                   '2048', '4', '.05', '1', '1.7', 'host', 'g5-series', 'no-neighbors', '64', 'candidate-only']
            log = out/f'{repeat}-{key}.log'; observations = []
            print('RUN', repeat, key, flush=True)
            with log.open('w') as f:
                child = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
                start = time.monotonic()
                while child.poll() is None:
                    if time.monotonic()-start > 300:
                        child.kill(); child.wait(); raise TimeoutError(key)
                    try:
                        threads = len(list(Path(f'/proc/{child.pid}/task').iterdir()))
                        freqs = {}
                        for c in cores:
                            try:
                                freqs[c] = int(Path(f'/sys/devices/system/cpu/cpu{c}/cpufreq/scaling_cur_freq').read_text())
                            except (FileNotFoundError, PermissionError):
                                freqs[c] = None
                        observations.append(dict(elapsed=time.monotonic()-start, threads=threads, khz=freqs))
                    except (FileNotFoundError, ProcessLookupError):
                        pass
                    time.sleep(.05)
            if child.returncode: raise RuntimeError(log.read_text()[-3000:])
            s = log.read_text()
            samples = [float(l.split()[3]) for l in s.splitlines() if l.startswith('TIMING ')]
            errors = [float(l.split()[4]) for l in s.splitlines() if l.startswith('CASE ')]
            if len(samples) != 5 or len(errors) != 1: raise RuntimeError('incomplete output')
            records[key]['samples'] += samples; records[key]['errors'] += errors
            run = dict(command=cmd, environment=controls, samples=samples, observations=observations,
                       median_seconds=statistics.median(samples),
                       peak_threads=max((o['threads'] for o in observations), default=None))
            records[key]['runs'].append(run)
            print('RESULT', repeat, key, run['median_seconds']*1000, 'peak_threads', run['peak_threads'], flush=True)
            (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    for record in records.values(): record['seconds'] = statistics.median(record['samples'])
    report['passed'] = True
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
