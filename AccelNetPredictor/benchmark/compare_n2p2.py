#!/usr/bin/env python3
"""Single-core n2p2/common-CPU and real-device GPU benchmark with result checks.

Use an OpenMP-disabled CPU build and a serial n2p2 reference. GPU subprocesses
inherit CUDA_VISIBLE_DEVICES and OMP_TARGET_OFFLOAD=MANDATORY. Report five warmed
wall-clock samples including forces (and AccelNet virial); initialization excluded.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import re
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'test'))
from check_n2p2_extended import run,compare
from n2p2_extended_fixtures import make_model,geometry,write_structure,TYPES


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('reference','cpu','gpu','output'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--atoms',type=int,nargs='+',default=[64,512])
    p.add_argument('--seconds',type=float,default=.1)
    p.add_argument('--scope',choices=('fixed','full'),default='fixed')
    p.add_argument('--affinity',type=int,default=6)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--max-slowdown',type=float,default=1.10)
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    if args.seconds <= 0:raise ValueError('--seconds must be positive')
    runtime_audit={}
    for key in ('reference','cpu'):
        symbols=subprocess.check_output(['nm','-u',str(getattr(args,key))],text=True)
        hits=re.findall(r'\b(?:GOMP_|__kmpc_|__nvomp|__pgi_omp|omp_get_|omp_set_)\w*',symbols)
        if hits:raise ValueError(f'{key} imports OpenMP runtime symbols: {hits[:5]}')
        runtime_audit[key]=dict(openmp_runtime_symbols=hits)
    os.sched_setaffinity(0,{args.affinity})
    os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY')
    rows=[];commands=[]
    def timed(cmd,tag):
        commands.append(list(map(str,cmd)));r=run(cmd,out/(tag+'.log'))
        samples=[x[0] for x in r['TIMING']]
        return r,dict(seconds=statistics.median(samples),samples=samples)
    for n in args.atoms:
        pos,cell=geometry(n);stem=out/f'structure-{n}';write_structure(stem,pos,cell)
        for kind in (9,*TYPES):
            model=out/f'type{kind}';make_model(model,(kind,),count=3 if kind==9 else 6)
            tag=f'{n}-type{kind}'
            ref,rt=timed([args.reference,model,stem.with_suffix('.data'),args.seconds,args.scope],tag+'-reference')
            cpu,ct=timed([args.cpu,model,stem.with_suffix('.xsf'),args.seconds,'cpu',args.scope,'n2p2'],tag+'-cpu')
            gpu,gt=timed([args.gpu,model,stem.with_suffix('.xsf'),args.seconds,'gpu',args.scope,'n2p2'],tag+'-gpu')
            row=dict(atoms=n,type=kind,n2p2=rt,cpu=ct,gpu=gt,errors=dict(cpu=compare(cpu,ref),gpu=compare(gpu,ref)),
                     cpu_vs_n2p2=rt['seconds']/ct['seconds'],gpu_vs_n2p2=rt['seconds']/gt['seconds'],gpu_vs_cpu=ct['seconds']/gt['seconds'])
            if kind==9:
                row['g5_modes']={}
                for backend,exe in [('host',args.cpu),('gpu',args.gpu)]:
                    for mode in (1,3):
                        r,t=timed([exe,model,stem.with_suffix('.xsf'),args.seconds,backend,args.scope,'n2p2',mode],tag+f'-{backend}-mode{mode}')
                        compare(r,ref);row['g5_modes'][f'{backend}-{mode}']=t
            rows.append(row);print(tag,row['cpu_vs_n2p2'],row['gpu_vs_n2p2'],flush=True)
            report=dict(scope=args.scope,affinity=args.affinity,runtime_audit=runtime_audit,rows=rows,commands=commands,
                        gpu_uuid=os.environ.get('CUDA_VISIBLE_DEVICES'),
                        binaries={key:dict(path=str(getattr(args,key)),sha256=hashlib.sha256(getattr(args,key).read_bytes()).hexdigest()) for key in ('reference','cpu','gpu')})
            (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    if args.baseline:
        baseline=json.loads(args.baseline.read_text())
        if baseline['scope']!=args.scope:raise ValueError('baseline scope differs')
        if baseline.get('gpu_uuid')!=report['gpu_uuid']:raise ValueError('baseline GPU differs')
        old={(r['atoms'],r['type']):r for r in baseline['rows']}
        for row in rows:
            before=old[(row['atoms'],row['type'])]
            for backend in ('cpu','gpu'):
                if row[backend]['seconds']>args.max_slowdown*before[backend]['seconds']:
                    raise AssertionError(f"performance regression: {row['atoms']} type {row['type']} {backend}")

if __name__=='__main__':main()
