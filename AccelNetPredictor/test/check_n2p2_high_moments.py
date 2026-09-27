#!/usr/bin/env python3
"""Check exact high-order G5 moments against unmodified n2p2, including mixtures."""
import argparse
import json
from pathlib import Path
from check_n2p2_extended import run, compare
from n2p2_extended_fixtures import TYPES, make_model, geometry, write_structure


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('reference', 'candidate', 'converter', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--backend', choices=('cpu', 'host', 'gpu'), default='cpu')
    a = p.parse_args()
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    degrees = (1, 10, 11, 12, 13, 14, 15, 16, 17, 11.5, 11 + 4e-13, 2.3)
    results = []
    for tag, kinds, elements in [('g5-high', (9,), ('O', 'H')),
                                  ('mixed-four', (2, 3, 9, *TYPES), ('O', 'N', 'C', 'H'))]:
        model = out / tag
        make_model(model, kinds, 'p2a', count=len(degrees), normalized=True, elements=elements)
        lines = []; counts = {}
        for line in (model / 'input.nn').read_text().splitlines():
            f = line.split()
            if f and f[0] == 'symfunction_short' and f[2] == '9':
                index = counts.get(f[1], 0); counts[f[1]] = index + 1
                f[7] = format(degrees[index % len(degrees)], '.17g')
                line = ' '.join(f)
            lines.append(line)
        (model / 'input.nn').write_text('\n'.join(lines) + '\n')
        positions, cell = geometry(16)
        cell = [[2*x for x in row] for row in cell]
        names = sorted(elements, key=lambda e: {'H': 1, 'C': 6, 'N': 7, 'O': 8}[e])
        stem = out / (tag + '-structure')
        write_structure(stem, positions, cell, [names[i % len(names)] for i in range(16)])
        ref = run([a.reference, model, stem.with_suffix('.data'), 0, 'fixed'], out / (tag + '-reference.log'))
        def evaluate(mode, native=False):
            return run([a.candidate, model / 'native' if native else model, stem.with_suffix('.xsf'),
                        0, a.backend, 'fixed', 'native' if native else 'n2p2', mode],
                       out / f'{tag}-mode{mode}{"-native" if native else ""}.log')
        values = {}
        for mode in (0, 1, 2, 3):
            values[mode] = evaluate(mode)
            errors = compare(values[mode], ref)
            if mode > 1: errors['direct'] = compare(values[mode], values[1])
            results.append(dict(case=tag, mode=mode, errors=errors))
        run([a.converter, 'n2p2-to-accelnet', model, model / 'native'], out / (tag + '-convert.log'))
        results.append(dict(case=tag+'-native', mode=3, errors=compare(evaluate(3, True), ref)))
        # Independent energy finite differences on the high-order mixture.
        h = 2e-5; force_error = 0.0; strain_error = 0.0
        for kind in ('force', 'strain'):
            for axis in range(3):
                energies = []
                for sign in (-1, 1):
                    pos = [v[:] for v in positions]; lat = [v[:] for v in cell]
                    if kind == 'force':
                        pos[0][axis] += sign*h
                    else:
                        for v in pos + lat: v[axis] += sign*h*v[(axis+1) % 3]
                    fd = out / f'{tag}-{kind}-{axis}-{sign}'
                    write_structure(fd, pos, lat, [names[i % len(names)] for i in range(16)])
                    r = run([a.reference, model, fd.with_suffix('.data'), 0, 'fixed'], fd.with_suffix('.log'))
                    energies.append(r['ENERGY'][0][0])
                expected = -(energies[1]-energies[0])/(2*h)
                actual = values[3]['FORCE'][0][axis] if kind == 'force' else values[3]['VIRIAL'][(axis+1) % 3][axis]
                error = abs(actual-expected)
                if error > 2e-8 + 2e-6*abs(expected): raise AssertionError((kind, error))
                if kind == 'force': force_error = max(force_error, error)
                else: strain_error = max(strain_error, error)
        results.append(dict(case=tag+'-finite-differences', force_error=force_error, strain_error=strain_error))
        print(tag, 'PASS', flush=True)
    report = dict(passed=True, backend=a.backend, degrees=degrees, cases=results)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
