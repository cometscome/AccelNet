import os, subprocess, pathlib, statistics, json
root=pathlib.Path('/tmp/accelnet-gpu-stage4')
env=os.environ.copy(); env.update(CUDA_VISIBLE_DEVICES='GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9',OMP_NUM_THREADS='1')
old=root/'before/accelnet-target-benchmark'
new=pathlib.Path('/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark')
results=[]
for n,p in [(8,3),(64,3),(512,3),(1024,3),(512,8)]:
    for label,binary in [('before',old),('after',new)]:
        command=['taskset','-c','6',str(binary),str(n),str(p),'0.2']
        output=subprocess.check_output(command,env=env,text=True,stderr=subprocess.STDOUT)
        (root/f'final-{label}-{n}-{p}.log').write_text(output)
        times=[list(map(float,line.split()[2:])) for line in output.splitlines() if line.startswith('TIMING ')]
        medians=[statistics.median(row[i] for row in times) for i in range(4)]
        case=next(line for line in output.splitlines() if line.startswith('CASE ')).split()
        results.append(dict(label=label,natoms=n,order=p,edges=int(case[3]),max_error=float(case[4]),median_seconds=medians))
        print(label,n,p,medians,flush=True)
for mode in [1,2,0]:
    output=subprocess.check_output(['taskset','-c','6',str(new),'512','8','0.2',str(mode),'1.1'],env=env,text=True,stderr=subprocess.STDOUT)
    (root/f'final-dense-mode{mode}.log').write_text(output)
    times=[list(map(float,line.split()[2:])) for line in output.splitlines() if line.startswith('TIMING ')]
    medians=[statistics.median(row[i] for row in times) for i in range(4)]
    results.append(dict(label='dense',mode=mode,natoms=512,order=8,spacing=1.1,median_seconds=medians))
    print('dense',mode,medians,flush=True)
(root/'final-performance.json').write_text(json.dumps(results,indent=2)+'\n')
