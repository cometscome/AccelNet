#!/usr/bin/env python3
"""Check unmodified trained n2p2 examples, including multi-element compact models.

Models remain external: --examples is n2p2/examples/nnp-predict. No downloads or
model redistribution are performed. Every consumed input is SHA256 recorded.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from check_n2p2_extended import run, compare

MODELS = ('Ethylbenzene_SCAN', 'Anisole_SCAN', 'DMABN_SCAN', 'H2O_RPBE-D3')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('reference','candidate','converter','examples','output'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--backend',choices=('cpu','host','gpu'),default='cpu')
    p.add_argument('--models',nargs='+',default=MODELS)
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    records=[]
    for name in args.models:
        source=args.examples.resolve()/name; dest=out/name;dest.mkdir(exist_ok=True)
        blocks=source.joinpath('input.data').read_text().split('end')
        sf=Counter()
        for line in source.joinpath('input.nn').read_text().splitlines():
            x=line.split()
            if x and x[0]=='symfunction_short':sf[(x[1],int(x[2]))]+=1
        native=dest/'native';exported=dest/'exported'
        run([args.converter,'n2p2-to-accelnet',source,native],dest/'convert.log')
        filenames=[native/line.split()[1] for line in (native/'networks.list').read_text().splitlines()]
        run([args.converter,'accelnet-to-n2p2',exported,*filenames],dest/'export.log')
        hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in source.iterdir()
                if f.name in ('input.nn','input.data','scaling.data') or f.name.startswith('weights.')}
        for frame,block in enumerate(blocks):
            if not any(line.strip()=='begin' for line in block.splitlines()):continue
            rows=[line.split() for line in block.splitlines()]
            atoms=[x for x in rows if x and x[0]=='atom'];lattice=[x[1:] for x in rows if x and x[0]=='lattice']
            stem=dest/f'frame-{frame}'
            stem.with_suffix('.data').write_text(block.strip()+'\nend\n')
            xsf=['CRYSTAL','PRIMVEC',*[' '.join(x) for x in lattice]] if lattice else []
            xsf+=['PRIMCOORD',f'{len(atoms)} 1',*[' '.join([x[4],*x[1:4]]) for x in atoms]]
            stem.with_suffix('.xsf').write_text('\n'.join(xsf)+'\n')
            ref=run([args.reference,source,stem.with_suffix('.data'),0,'fixed'],dest/f'{frame}-reference.log')
            errors={}
            for label,directory,fmt in [('import',source,'n2p2'),('native',native,'native')]:
                new=run([args.candidate,directory,stem.with_suffix('.xsf'),0,args.backend,'fixed',fmt],dest/f'{frame}-{label}.log')
                errors[label]=compare(new,ref)
            new=run([args.reference,exported,stem.with_suffix('.data'),0,'fixed'],dest/f'{frame}-roundtrip.log')
            errors['roundtrip']=compare(new,ref)
            record=dict(model=name,frame=frame,atoms=len(atoms),elements=sorted(set(x[4] for x in atoms)),
                        descriptors={f'{el}:{kind}':count for (el,kind),count in sorted(sf.items())},errors=errors,sha256=hashes)
            records.append(record);print(name,frame,'PASS',errors,flush=True)
    (out/'report.json').write_text(json.dumps(dict(passed=True,backend=args.backend,cases=records),indent=2)+'\n')

if __name__=='__main__':main()
