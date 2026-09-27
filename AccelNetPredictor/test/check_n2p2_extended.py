#!/usr/bin/env python3
"""Independent n2p2 checks for all extended SF types, conversion and strain derivatives."""
import argparse
import json
import math
from pathlib import Path
import subprocess
from n2p2_extended_fixtures import TYPES, SUBTYPES, make_model, geometry, write_structure


def run(command, log):
    p = subprocess.run(list(map(str,command)), capture_output=True, text=True, timeout=180)
    log.write_text(p.stdout+p.stderr)
    if p.returncode:
        raise RuntimeError(f'{command}: {p.stderr[-1500:]}')
    result={}
    for line in p.stdout.splitlines():
        fields=line.split()
        if fields and fields[0] in ('ENERGY','FORCE','VIRIAL','TIMING'):
            result.setdefault(fields[0],[]).append([float(x) for x in fields[1:]])
    return result


def compare(a,b):
    errors={}
    for key in ('ENERGY','FORCE'):
        av=[v for row in a[key] for v in row];bv=[v for row in b[key] for v in row]
        if len(av)!=len(bv): raise AssertionError('shape mismatch')
        errors[key]=max(abs(x-y) for x,y in zip(av,bv))
        if any(not math.isfinite(x) or not math.isfinite(y) or abs(x-y)>2e-10+2e-9*abs(y) for x,y in zip(av,bv)):
            raise AssertionError(f'{key}: {errors[key]}')
    return errors


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference',type=Path)
    ap.add_argument('--candidate',type=Path,required=True)
    ap.add_argument('--converter',type=Path,required=True)
    ap.add_argument('--backend',choices=('cpu','host','gpu'),default='cpu')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    args.reference=args.reference.resolve() if args.reference else None;args.candidate=args.candidate.resolve();args.converter=args.converter.resolve()
    positions,cell=geometry(8)
    # Keep cell vectors longer than Rc: compact windows centered on 0/pi have
    # an upstream n2p2 discontinuity at collinear periodic self-image pairs.
    # Exercise exact axial collinearity separately below.
    cell=[[1.5*x for x in row] for row in cell]
    write_structure(out/'structure',positions,cell)
    cases=[(f'type{k}-{sub}',(k,),sub,1,False) for k in TYPES for sub in (('p2',) if k<20 else SUBTYPES)]
    cases += [(f'weighted-cutoff{ct}',(12,13),'p2',ct,False) for ct in range(9)]
    cases += [('mixed-normalized',(2,3,9,*TYPES),'p3a',7,True)]
    def reference(directory,stem,log):
        if args.reference:
            return run([args.reference,directory,stem.with_suffix('.data'),0,'fixed'],log)
        return run([args.candidate,directory,stem.with_suffix('.xsf'),0,'cpu','fixed','n2p2'],log)
    results=[]
    for tag,types,sub,ct,norm in cases:
        directory=out/tag;make_model(directory,types,sub,cutoff=ct,normalized=norm)
        ref=reference(directory,out/'structure',out/(tag+'-reference.log'))
        new=run([args.candidate,directory,out/'structure.xsf',0,args.backend,'fixed','n2p2'],out/(tag+'-candidate.log'))
        errors=compare(new,ref)
        native=directory/'native'
        run([args.converter,'n2p2-to-accelnet',directory,native],out/(tag+'-convert.log'))
        converted=run([args.candidate,native,out/'structure.xsf',0,args.backend,'fixed','native'],out/(tag+'-native.log'))
        errors['converted']=compare(converted,ref)
        exported=directory/'exported'
        run([args.converter,'accelnet-to-n2p2',exported,native/'H.nn.ascii',native/'O.nn.ascii'],out/(tag+'-export.log'))
        roundtrip=reference(exported,out/'structure',out/(tag+'-roundtrip.log'))
        errors['roundtrip']=compare(roundtrip,ref)
        results.append(dict(case=tag,errors=errors))
        print(tag,'PASS',errors,flush=True)
    # Independent strain finite differences on a mixture of every supported type.
    directory=out/'mixed-normalized';h=2e-5
    base=run([args.candidate,directory,out/'structure.xsf',0,args.backend,'fixed','n2p2'],out/'strain-base.log')
    max_strain=0
    for a in range(3):
        for b in range(3):
            energies=[]
            for sign in (-1,1):
                transform=lambda vectors:[[v[k]+(sign*h*v[b] if k==a else 0) for k in range(3)] for v in vectors]
                stem=out/f'strain-{a}-{b}-{sign}'
                write_structure(stem,transform(positions),transform(cell))
                r=reference(directory,stem,stem.with_suffix('.log'))
                energies.append(r['ENERGY'][0][0])
            expected=-(energies[1]-energies[0])/(2*h)
            actual=base['VIRIAL'][b][a]
            max_strain=max(max_strain,abs(actual-expected))
            if abs(actual-expected)>2e-8+2e-6*abs(expected): raise AssertionError(f'strain {a} {b}: {actual} {expected}')
    # Isolated atoms, nonperiodic and exact collinear inputs exercise sparse/endpoint cases.
    for tag,pos in [('isolated',[[0,0,0]]),('collinear',[[0,0,0],[1,0,0],[-1.2,0,0]]),('cluster',positions)]:
        stem=out/tag;write_structure(stem,pos,None)
        ref=reference(directory,stem,out/(tag+'-reference.log'))
        new=run([args.candidate,directory,stem.with_suffix('.xsf'),0,args.backend,'fixed','n2p2'],out/(tag+'-candidate.log'))
        results.append(dict(case=tag,errors=compare(new,ref)))
    # Finite difference of the independent n2p2 energy versus all shared force components.
    force_error=0
    for atom in range(len(positions)):
        for axis in range(3):
            energies=[]
            for sign in (-1,1):
                moved=[v[:] for v in positions];moved[atom][axis]+=sign*h
                stem=out/f'force-{atom}-{axis}-{sign}';write_structure(stem,moved,cell)
                r=reference(directory,stem,stem.with_suffix('.log'))
                energies.append(r['ENERGY'][0][0])
            expected=-(energies[1]-energies[0])/(2*h);actual=base['FORCE'][atom][axis]
            force_error=max(force_error,abs(actual-expected))
            if abs(actual-expected)>2e-8+2e-6*abs(expected): raise AssertionError('force finite difference')
    # Four species, all eleven SF types, ten channels/type, including every
    # unordered chemical pair. Reverse-declared elements test canonical sorting.
    four=out/'four-elements';make_model(four,(2,3,9,*TYPES),'p2a',count=10,normalized=True,elements=('O','N','C','H'))
    pos,lat=geometry(16);lat=[[2*x for x in row] for row in lat]
    stem=out/'four-structure';write_structure(stem,pos,lat,['H','C','N','O']*4)
    ref=reference(four,stem,out/'four-reference.log')
    for mode in (0,1,3):
        backend=args.backend if mode==0 or args.backend=='gpu' else 'host'
        new=run([args.candidate,four,stem.with_suffix('.xsf'),0,backend,'fixed','n2p2',mode],out/f'four-mode{mode}.log')
        results.append(dict(case=f'four-elements-mode{mode}',errors=compare(new,ref)))
    report=dict(passed=True,independent_n2p2=bool(args.reference),backend=args.backend,cases=results,max_strain_error=max_strain,max_force_fd_error=force_error)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',len(results),'cases; strain error',max_strain,'force FD',force_error)

if __name__=='__main__':main()
