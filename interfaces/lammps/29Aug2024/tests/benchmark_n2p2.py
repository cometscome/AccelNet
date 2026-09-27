#!/usr/bin/env python3
"""End-to-end LAMMPS CPU/GPU comparison on unmodified n2p2 prediction models.

One process/core, OpenMP disabled in the CPU binaries. Compare energies, forces,
virial pressure and short NVE trajectories before accepting timings. All timed
segments advance the same trajectory. Loading, run setup and dumps are excluded
from LAMMPS Loop time; integration and in-loop neighbor rebuilds are included.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import time
import numpy as np

MODELS=('Ethylbenzene_SCAN','Anisole_SCAN','DMABN_SCAN','H2O_RPBE-D3')
MASS={'H':1.00794,'C':12.0107,'N':14.0067,'O':15.9994}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p):
    lines=p.read_text().splitlines()
    start=next(i for i,l in enumerate(lines) if l.startswith('ITEM: ATOMS'))+1
    return np.array([list(map(float,l.split())) for l in lines[start:]])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('before','after','converter','examples','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--gpu',type=Path)
    p.add_argument('--gpu-before',type=Path)
    p.add_argument('--models',nargs='+',default=MODELS)
    p.add_argument('--samples',type=int,default=3)
    p.add_argument('--rounds',type=int,default=2)
    p.add_argument('--variants',nargs='+',choices=('n2p2','before','after','gpu_before','gpu'),
                   default=['n2p2','before','after','gpu_before','gpu'])
    p.add_argument('--steps',type=int,help='Default: 100 for isolated DMABN, 10 for water, 3 for dense compact models')
    p.add_argument('--affinity',type=int,default=6)
    p.add_argument('--max-n2p2-slowdown',type=float,
                   help='Reject an after/n2p2 CPU step-time ratio above this limit')
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    os.sched_setaffinity(0,{args.affinity})
    # Warm the pinned core out of its idle frequency state before comparisons.
    warmup_end = time.monotonic() + 20.0
    while time.monotonic() < warmup_end:
        sum(range(10000))
    os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY')
    binaries={}
    for name in ('before','after','gpu','gpu_before'):
        exe=getattr(args,name)
        if exe is None:continue
        exe=exe.resolve();setattr(args,name,exe)
        if name in ('before','after'):
            symbols=subprocess.check_output(['nm','-u',exe],text=True)
            if re.search(r'\b(?:GOMP_|__kmpc_|__nvomp|__pgi_omp|omp_get_|omp_set_)\w*',symbols):
                raise ValueError('CPU binary imports OpenMP: '+str(exe))
        binaries[name]=dict(path=str(exe),sha256=sha(exe))
    report=dict(binaries=binaries,affinity=args.affinity,samples=args.samples,rounds=args.rounds,units='electron',
                timestep_fs=.1,neighbor='skin .6; every 5 delay 0 check no',
                gpu_uuid=os.environ.get('CUDA_VISIBLE_DEVICES'),core_warmup_seconds=20,rows=[])
    for name in args.models:
        model=args.examples.resolve()/name
        rows=[l.split() for l in (model/'input.data').read_text().splitlines()]
        atoms=[r for r in rows if r and r[0]=='atom']
        lattice=np.array([list(map(float,r[1:])) for r in rows if r and r[0]=='lattice'])
        case=out/name;case.mkdir(exist_ok=True)
        native=case/'native'
        converted=subprocess.run([args.converter.resolve(),'n2p2-to-accelnet',model,native],capture_output=True,text=True,check=True)
        (case/'conversion.log').write_text(converted.stdout+converted.stderr)
        elements=[line.split()[0] for line in (native/'networks.list').read_text().splitlines()]
        xyz=np.array([list(map(float,r[1:4])) for r in atoms])
        settings=[l.split('#')[0].split() for l in (model/'input.nn').read_text().splitlines()]
        sf=[r for r in settings if r and r[0]=='symfunction_short']
        rc=max(float(r[{2:6,3:8,20:5,22:6}[int(r[2])]]) for r in sf)
        periodic=lattice.size>0
        if periodic:
            if not np.allclose(lattice,np.diag(np.diag(lattice))):raise ValueError('orthogonal examples required')
            lo=np.zeros(3);hi=np.diag(lattice)
        else:
            lo=xyz.min(axis=0)-rc-2;hi=xyz.max(axis=0)+rc+2
        data=f'n2p2 official {name}\n\n{len(atoms)} atoms\n{len(elements)} atom types\n\n'
        for d,axis in enumerate('xyz'): data+=f'{lo[d]:.17g} {hi[d]:.17g} {axis}lo {axis}hi\n'
        data+='\nMasses\n\n'+''.join(f'{i+1} {MASS[e]}\n' for i,e in enumerate(elements))
        data+='\nAtoms # atomic\n\n'+''.join(f'{i+1} {elements.index(at[4])+1} '+' '.join(at[1:4])+'\n' for i,at in enumerate(atoms))
        (case/'structure.data').write_text(data)
        steps=args.steps or (100 if not periodic else 10 if len(atoms)>1000 else 3)
        paths=[('n2p2',args.before),('before',args.before),('after',args.after)]
        if args.gpu_before:paths.append(('gpu_before',args.gpu_before))
        if args.gpu:paths.append(('gpu',args.gpu))
        paths=[item for item in paths if item[0] in args.variants]
        if not paths or paths[0][0]!='n2p2':raise ValueError('n2p2 must be included as the independent reference')
        records={}; reference=None
        for repeat in range(args.rounds):
            ordered=paths if repeat%2==0 else list(reversed(paths))
            for label,exe in ordered:
                work=case/(label+f'-{repeat}');work.mkdir(exist_ok=True)
                gpu=label.startswith('gpu')
                package='package gpu 1 neigh yes newton on split 1' if gpu else ''
                if label=='n2p2':
                    pair=f'pair_style hdnnp {rc+1e-10:.17g} dir {model} showew no resetew no showewsum 0 maxew -1 cflength 1 cfenergy 1\npair_coeff * * '+' '.join(elements)
                elif gpu:
                    pair='pair_style accelnet/gpu auto '+' '.join(str(native/(e+'.nn.ascii')) for e in elements)+'\npair_coeff * *'
                else:
                    pair=f'pair_style accelnet auto n2p2 {model} '+' '.join(elements)+'\npair_coeff * *'
                text=f'''clear
    {package}
    units electron
    atom_style atomic
    boundary {'p p p' if periodic else 'f f f'}
    read_data {case}/structure.data
    {pair}
    neighbor 0.6 bin
    neigh_modify every 5 delay 0 check no one 10000 page 1000000
    atom_modify sort 0 0.0
    velocity all create 100 4928459 mom yes rot no dist gaussian
    fix integrate all nve
    timestep 0.1
    compute vir all pressure NULL virial
    thermo {steps}
    thermo_style custom step pe c_vir[1] c_vir[2] c_vir[3] c_vir[4] c_vir[5] c_vir[6]
    thermo_modify format float %.17g lost error
    run 0
    write_dump all custom initial.dump id type x y z fx fy fz modify sort id format float %.17g
    run 2
    '''+f'run {steps}\n'*args.samples+'''write_dump all custom final.dump id type x y z fx fy fz modify sort id format float %.17g
    '''
                (work/'in.benchmark').write_text(text)
                command=[str(exe),'-in','in.benchmark']
                print('RUN',name,label,flush=True)
                with (work/'stdout.txt').open('w') as f:
                    proc=subprocess.run(command,cwd=work,stdout=f,stderr=subprocess.STDOUT,timeout=900)
                log=(work/'stdout.txt').read_text()
                if proc.returncode:raise RuntimeError(str(work)+'\n'+log[-2000:])
                loops=[float(v)/steps for v in re.findall(r'Loop time of ([\d.eE+-]+)',log)[-args.samples:]]
                pairs=[float(v)/steps for v in re.findall(r'^Pair\s+\|\s*[\d.eE+-]+\s*\|\s*([\d.eE+-]+)',log,re.M)[-args.samples:]]
                thermo=[]
                for line in log.splitlines():
                    f=line.split()
                    if len(f)==8 and re.fullmatch(r'\d+',f[0]):
                        try:thermo.append(list(map(float,f)))
                        except ValueError:pass
                values=dict(initial=dump(work/'initial.dump'),final=dump(work/'final.dump'),thermo=np.array(thermo))
                if reference is None:reference=values
                errors={}
                for key in ('initial','final'):
                    np.testing.assert_allclose(values[key][:,:2],reference[key][:,:2],rtol=0,atol=0)
                    np.testing.assert_allclose(values[key][:,2:],reference[key][:,2:],rtol=2e-9,atol=2e-8)
                    errors[key+'_force']=float(np.max(np.abs(values[key][:,5:]-reference[key][:,5:])))
                    errors[key+'_position']=float(np.max(np.abs(values[key][:,2:5]-reference[key][:,2:5])))
                np.testing.assert_allclose(values['thermo'][:,:2],reference['thermo'][:,:2],rtol=2e-9,atol=2e-8)
                np.testing.assert_allclose(values['thermo'][:,2:],reference['thermo'][:,2:],rtol=2e-8,atol=.05)
                errors['energy']=float(np.max(np.abs(values['thermo'][:,1]-reference['thermo'][:,1])))
                errors['virial_pressure_pa']=float(np.max(np.abs(values['thermo'][:,2:]-reference['thermo'][:,2:])))
                if len(loops)!=args.samples or len(pairs)!=args.samples:raise ValueError('missing timing samples')
                if label not in records:records[label]=dict(samples=[],pair_samples=[],errors=[],command=command)
                records[label]['samples'].extend(loops);records[label]['pair_samples'].extend(pairs)
                records[label]['errors'].append(errors)
                records[label]['seconds_per_step']=statistics.median(records[label]['samples'])
                print('RESULT',name,label,records[label]['seconds_per_step'],errors,flush=True)
        record=dict(model=name,atoms=len(atoms),elements=elements,steps=steps,records=records,
                    model_sha256={f.name:sha(f) for f in model.iterdir() if f.name in ('input.nn','input.data','scaling.data') or f.name.startswith('weights.')})
        report['rows'].append(record);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    if args.max_n2p2_slowdown is not None:
        if 'after' not in args.variants:
            raise ValueError('--max-n2p2-slowdown requires the after variant')
        report['max_n2p2_slowdown']=args.max_n2p2_slowdown
        report['n2p2_regressions']=[
            (row['model'],row['records']['after']['seconds_per_step']/row['records']['n2p2']['seconds_per_step'])
            for row in report['rows']
            if row['records']['after']['seconds_per_step'] >
                args.max_n2p2_slowdown*row['records']['n2p2']['seconds_per_step']
        ]
        if report['n2p2_regressions']:
            (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
            raise AssertionError(f"n2p2 parity check: {report['n2p2_regressions']}")
    report['passed']=True;(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
