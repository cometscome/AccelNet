from pathlib import Path
import subprocess,os,json
out=Path('/tmp/accelnet-g4-tuning')
os.sched_setaffinity(0,{6})
env=dict(os.environ,OMP_NUM_THREADS='1',OMP_TARGET_OFFLOAD='MANDATORY',CUDA_VISIBLE_DEVICES='GPU-2644154d-7268-af42-6631-59e1f3c6e7f3')
for tag in ['before','scalars','pair-parallel']:
    binary=out/tag/'gpu'
    cmd=['/usr/local/cuda-13.1/bin/ncu','--kernel-name','regex:g4_values_derivatives','--launch-count','1',
         '--section','LaunchStats','--section','Occupancy','--section','SpeedOfLight','--section','MemoryWorkloadAnalysis',
         '--export',str(out/f'profile-{tag}'),str(binary),'4096','8','.001','1','1.7','gpu','g4-series','no-neighbors']
    try:
        p=subprocess.run(cmd,env=env,text=True,capture_output=True,timeout=180)
        (out/f'profile-{tag}.log').write_text(p.stdout+p.stderr)
        (out/f'profile-{tag}-command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        print(tag,p.returncode,p.stdout[-2000:],p.stderr[-1000:],flush=True)
        if 'ERR_NVGPUCTRPERM' in p.stdout+p.stderr:break
    except subprocess.TimeoutExpired:
        print('TIMEOUT',tag,flush=True);break
