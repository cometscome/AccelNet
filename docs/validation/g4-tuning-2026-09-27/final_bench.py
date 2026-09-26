from pathlib import Path
import json,shutil
from run import out,bench,uuids
for label in ['h100','blackwell','gnu','nvhpc']:
    compiler='gpu' if label in uuids else label
    current=Path('/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark') if compiler=='gpu' else Path(f'/tmp/accelnet-other-descriptors/build-{compiler}-serial/bin/accelnet-target-benchmark')
    saved=out/'final';saved.mkdir(exist_ok=True);shutil.copy2(current,saved/compiler)
    variants=[('v13',Path('/tmp/accelnet-g4-fused/before')/compiler),('v14',out/'before'/compiler),('v15',saved/compiler)]
    rows=[]
    for n in [512,4096]:
      for family in ['g4-distinct','g4-series']:
        for rep in range(2):
          for variant,binary in variants[::1 if rep==0 else -1]:
            result=bench(f'final-{variant}-{rep}',label,binary,sizes=(n,),families=(family,),seconds='.10')
            for row in result:row['round']=rep;row['variant']=variant
            rows+=result;(out/f'final-{label}.json').write_text(json.dumps(rows,indent=2)+'\n')
