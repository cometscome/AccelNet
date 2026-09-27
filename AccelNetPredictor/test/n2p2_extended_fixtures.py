"""Deterministic synthetic n2p2 models; these are numerical fixtures, not trained potentials."""
from pathlib import Path
import math
import random

TYPES = (12, 13, 20, 21, 22, 23, 24, 25)
SUBTYPES = ('e', 'p1', 'p1a', 'p2', 'p2a', 'p3', 'p3a', 'p4', 'p4a')

def make_model(directory, types=TYPES, subtype='p2', count=6, cutoff=1, normalized=False, elements=('H','O')):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    numbers={'H':1,'C':6,'N':7,'O':8}
    elements=tuple(sorted(elements,key=numbers.__getitem__))
    pairs=[(a,b) for i,a in enumerate(elements) for b in elements[i:]]
    lines = [f'number_of_elements {len(elements)}', 'elements '+' '.join(reversed(elements)), 'nnp_type 2G-HDNNP',
             f'cutoff_type {cutoff} 0.15', 'scale_symmetry_functions',
             'scale_min_short -1', 'scale_max_short 1', 'global_hidden_layers_short 2',
             'global_nodes_short 8 4', 'global_activation_short t t l']
    lines += [f'atom_energy {element} {0.03-.01*i}' for i,element in enumerate(elements)]
    if normalized:
        lines += ['mean_energy 0.17', 'conv_energy 1.7', 'conv_length 1.3', 'normalize_nodes']
    for species in elements:
        functions = []
        for kind in types:
            for i in range(count):
                rc = 3.8 + .1 * (i % 2)
                rl = -.4 + .2 * (i % 3)
                eta = .08 + .06 * (i % 3)
                rs = .12 * (i % 2)
                lam = (-1, 1)[i % 2]
                zeta = (1, 2, 4, 1.5)[i % 4]
                a, b = pairs[i % len(pairs)]
                left, right = ((0,180), (20,160), (-60,60), (120,240))[i % 4]
                if kind == 12: params = f'{eta} {rs} {rc}'
                elif kind == 13: params = f'{eta} {rs} {lam} {zeta} {rc}'
                elif kind == 20: params = f'{a} {rl} {rc} {subtype}'
                elif kind in (21,22): params = f'{a} {b} {rl} {rc} {left} {right} {subtype}'
                elif kind == 23: params = f'{rl} {rc} {subtype}'
                elif kind in (24,25): params = f'{rl} {rc} {left} {right} {subtype}'
                elif kind == 2: params = f'{a} {eta} {rs} {rc}'
                elif kind in (3,9): params = f'{a} {b} {eta} {lam} {zeta} {rc} {rs}'
                else: raise ValueError(kind)
                functions.append(f'symfunction_short {species} {kind} {params}')
        random.Random(7281).shuffle(functions)
        lines.extend(functions)
        n = len(functions); nodes=(n,8,4,1)
        size = sum((a+1)*b for a,b in zip(nodes,nodes[1:]))
        weights = [.09*math.sin((i+1)*1.37 + (0 if species=='H' else .4)) for i in range(size)]
        (directory/f'weights.{numbers[species]:03d}.data').write_text(''.join(f'{w:.17g}\n' for w in weights))
    (directory/'input.nn').write_text('\n'.join(lines)+'\n')
    (directory/'scaling.data').write_text(''.join(
        f'{s} {i+1} 0 {100+.7*i:.17g} 1 10\n' for s in range(1,len(elements)+1) for i in range(n)))

def geometry(n, periodic=True):
    side=math.ceil(n**(1/3)-1e-12)
    positions=[]
    for i in range(n):
        positions.append([1.65*(i%side)+.031*math.sin(i*1.3),
                          1.65*((i//side)%side)+.027*math.cos(i*.9),
                          1.65*(i//(side*side))+.023*math.sin(i*.7)])
    cell=[[1.65*side,0,0],[.13,1.65*side,0],[.09,.17,1.65*side]] if periodic else None
    return positions,cell

def write_structure(stem, positions, cell, elements=None):
    stem=Path(stem);stem.parent.mkdir(parents=True,exist_ok=True)
    fmt=lambda row:' '.join(f'{x:.17g}' for x in row)
    xsf=[];data=['begin','comment n2p2 extension fixture']
    if cell is not None:
        xsf=['CRYSTAL','PRIMVEC',*[fmt(r) for r in cell]]
        data += ['lattice '+fmt(r) for r in cell]
    xsf+=['PRIMCOORD',f'{len(positions)} 1']
    for i,r in enumerate(positions):
        element=elements[i] if elements is not None else ('H','O')[i%2]
        xsf.append(element+' '+fmt(r));data.append('atom '+fmt(r)+' '+element+' 0 0 0 0 0')
    data+=['energy 0','charge 0','end']
    stem.with_suffix('.xsf').write_text('\n'.join(xsf)+'\n')
    stem.with_suffix('.data').write_text('\n'.join(data)+'\n')
