#!/usr/bin/env python3
"""Build an untouched n2p2 libnnp snapshot and benchmark, without OpenMP/MPI."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source',type=Path,required=True,help='n2p2 source root')
p.add_argument('--eigen',type=Path,required=True,help='directory containing Eigen/')
p.add_argument('--output',type=Path,required=True)
p.add_argument('--cxx',default='g++')
p.add_argument('--jobs',type=int,default=4)
a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
src=a.source.resolve()/'src/libnnp'
flags=['-O3','-std=c++11','-DEIGEN_DONT_PARALLELIZE','-I'+str(a.eigen.resolve()),'-I'+str(src)]
commands=[[a.cxx,*flags,'-c',str(f),'-o',str(out/(f.stem+'.o'))] for f in sorted(src.glob('*.cpp'))]
if not commands:raise ValueError('no libnnp sources')
with ThreadPoolExecutor(max_workers=a.jobs) as pool:
    list(pool.map(lambda c:subprocess.run(c,check=True),commands))
archive=['ar','rcs',str(out/'libnnp.a'),*[c[-1] for c in commands]]
subprocess.run(archive,check=True)
driver=Path(__file__).resolve().with_name('n2p2_reference_benchmark.cpp')
link=[a.cxx,*flags,str(driver),str(out/'libnnp.a'),'-o',str(out/'n2p2-reference')]
subprocess.run(link,check=True)
(out/'build.json').write_text(json.dumps(dict(commands=[*commands,archive,link],
    source_sha256={str(f.relative_to(a.source.resolve())):hashlib.sha256(f.read_bytes()).hexdigest()
                   for f in sorted(src.iterdir()) if f.suffix in ('.cpp','.h')},
    compiler=subprocess.check_output([a.cxx,'--version'],text=True)),indent=2)+'\n')
