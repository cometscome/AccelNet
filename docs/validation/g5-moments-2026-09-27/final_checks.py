from pathlib import Path
import subprocess,os,json,statistics
out=Path('/tmp/accelnet-g5-moment'); root=Path('/home/nagai/AccelNetGPU/AccelNet-clone')
os.sched_setaffinity(0,{6})
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3')
for compiler in ['gnu','nvhpc']:
 build=Path(f'/tmp/accelnet-other-descriptors/build-{compiler}-serial')
 for tag,cmd in [('configure',['cmake','-S',str(root),'-B',str(build),'-DACCELNET_TEST_G5_MOMENT_PERFORMANCE=ON']),('gate',['ctest','--test-dir',str(build),'-R','^predictor_g5_moment_performance$','--output-on-failure','-V'])]:
  with (out/f'moment-{tag}-{compiler}.log').open('w') as f:subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 (out/f'moment-gate-{compiler}.json').write_text((build/'AccelNetPredictor/g5-moment-performance.json').read_text())
 print('moment gate',compiler,'passed',flush=True)
rows=[]
for backend,exe,kind in [('gnu','/tmp/accelnet-other-descriptors/build-gnu-serial/bin/accelnet-target-benchmark','cpu-shared'),('nvhpc','/tmp/accelnet-other-descriptors/build-nvhpc-serial/bin/accelnet-target-benchmark','cpu-shared'),('h100','/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark','gpu')]:
 for rnd in range(2):
  for mode in ([1,3] if rnd==0 else [3,1]):
   cmd=[exe,'512','8','.05',str(mode),'1.7',kind,'g5-series','no-neighbors']
   p=subprocess.run(cmd,env=env,text=True,capture_output=True,check=True)
   tag=f'control-order8-{backend}-mode{mode}-round{rnd}';(out/(tag+'.log')).write_text(p.stdout+p.stderr)
   samples=[list(map(float,l.split()[2:4])) for l in p.stdout.splitlines() if l.startswith('TIMING ')]
   rows.append(dict(tag=tag,command=cmd,samples=samples,old_ms=1000*statistics.median(x[0] for x in samples),new_ms=1000*statistics.median(x[1] for x in samples)))
   (out/'order8-control.json').write_text(json.dumps(rows,indent=2)+'\n')
   print(tag,rows[-1]['old_ms'],rows[-1]['new_ms'],flush=True)
rows=[]
for compiler in ['gnu','nvhpc']:
 for family in ['lj','g4-distinct','g5-series']:
  for rnd in range(2):
   for variant in (['before','after'] if rnd==0 else ['after','before']):
    exe=str(out/f'before-{compiler}') if variant=='before' else f'/tmp/accelnet-other-descriptors/build-{compiler}-serial/bin/accelnet-target-benchmark'
    cmd=[exe,'512','8','.05','1','1.7','cpu-shared',family,'no-neighbors']
    p=subprocess.run(cmd,env=env,text=True,capture_output=True,check=True)
    tag=f'control-direct-{compiler}-{family}-{variant}-round{rnd}';(out/(tag+'.log')).write_text(p.stdout+p.stderr)
    samples=[list(map(float,l.split()[2:4])) for l in p.stdout.splitlines() if l.startswith('TIMING ')]
    rows.append(dict(compiler=compiler,family=family,variant=variant,round=rnd,command=cmd,samples=samples,new_ms=1000*statistics.median(x[1] for x in samples)))
    (out/'direct-control.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(tag,rows[-1]['new_ms'],flush=True)
