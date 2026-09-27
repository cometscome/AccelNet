from pathlib import Path
import os,sys,subprocess,json,statistics,hashlib
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-g4-tuning')
gpu='/tmp/accelnet-gpu-research/build-nvhpc'
uuids={'h100':'GPU-2644154d-7268-af42-6631-59e1f3c6e7f3','blackwell':'GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9'}
def build(label):
    path=gpu if label=='gpu' else f'/tmp/accelnet-other-descriptors/build-{label}-serial'
    with (out/f'build-{label}.log').open('w') as f:
        subprocess.run(['cmake','--build',path,'-j','8','--target','accelnet-target-benchmark','test_batch_target'],stdout=f,stderr=subprocess.STDOUT,check=True)
    return Path(path)/'bin/accelnet-target-benchmark'
def bench(tag,label,binary,sizes=(512,4096),families=('g4-distinct','g4-series'),seconds='.03'):
    os.sched_setaffinity(0,{6});rows=[]
    for n in sizes:
      for family in families:
        cmd=[str(binary),str(n),'8',seconds,'1','1.7','gpu' if label in uuids else 'host',family,'no-neighbors']
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES=uuids.get(label,''))
        p=subprocess.run(cmd,env=env,text=True,capture_output=True,timeout=600)
        (out/f'{tag}-{label}-{n}-{family}.log').write_text(p.stdout+p.stderr)
        if p.returncode:raise RuntimeError(p.stdout+p.stderr)
        samples=[list(map(float,l.split()[2:])) for l in p.stdout.splitlines() if l.startswith('TIMING')]
        phases=[list(map(float,l.split()[3:])) for l in p.stdout.splitlines() if l.startswith('PROFILE')]
        case=next(l.split() for l in p.stdout.splitlines() if l.startswith('CASE '))
        row=dict(variant=tag,backend=label,n=n,family=family,command=cmd,sha256=hashlib.sha256(Path(binary).read_bytes()).hexdigest(),samples=samples,phases=phases,max_error=float(case[4]),edges=int(case[3]),ms=1000*statistics.median(x[1] for x in samples),reference_ms=1000*statistics.median(x[0] for x in samples),ratio=statistics.median(x[1]/x[0] for x in samples))
        rows.append(row);print(tag,label,n,family,round(row['ms'],3),round(row['ratio'],3),flush=True)
        (out/f'{tag}-{label}.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows
if __name__=='__main__':
  if sys.argv[1]=='threads':
    import shutil
    source=root/'AccelNetPredictor/src/accelnet_target_kernels.f90';baseline=(out/'before/kernels.f90').read_text()
    for label in uuids:bench('before',label,out/'before/gpu')
    for limit in [32,64,128]:
      begin=baseline.index('    subroutine g4_values_derivatives(');i=baseline.index('if(device /= omp_get_initial_device()) &',begin)
      source.write_text(baseline[:i]+baseline[i:].replace('if(device /= omp_get_initial_device()) &',f'if(device /= omp_get_initial_device()) thread_limit({limit}) &',1))
      binary=build('gpu');saved=out/f'threads{limit}';shutil.copy2(binary,saved)
      for label in uuids:bench(f'threads{limit}',label,saved)
  else:
    tag,label=sys.argv[1:3];binary=build('gpu' if label in uuids else label)
    bench(tag,label,binary)
