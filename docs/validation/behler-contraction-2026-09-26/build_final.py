from pathlib import Path
import subprocess
out=Path('/tmp/accelnet-angular-contraction');root=Path('/home/nagai/AccelNetGPU/AccelNet-clone')
for name,build in [('gnu','/tmp/accelnet-other-descriptors/build-gnu-serial'),('nvhpc','/tmp/accelnet-other-descriptors/build-nvhpc-serial'),('gpu','/tmp/accelnet-gpu-research/build-nvhpc')]:
 cmd=['cmake','--build',build,'-j','8','--target','accelnet-target-benchmark','test_batch_target','test_batch','test_target_c','write_target_fixtures']
 with (out/f'build-final-{name}.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:raise SystemExit('FAIL '+name)
 print('BUILT',name,flush=True)
 baseline=out/('baseline-'+name);mods=baseline/'driver-modules';mods.mkdir(exist_ok=True)
 cc=['gfortran','-O3'] if name=='gnu' else ['/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers/bin/nvfortran','-fast','-O3']
 includes=[f'-I{baseline / d}' for d in ['AccelNetPredictor/modules','AccelNetPredictor/target-modules','AccelNetDescriptors/modules']]
 cc+=(['-J'+str(mods)] if name=='gnu' else ['-module',str(mods)])+includes
 if name=='gpu':cc+=['-mp=gpu','-gpu=cc90,cc120']
 cmd=cc+[str(root/'AccelNetPredictor/test/batch_test_support.f90'),str(root/'AccelNetPredictor/benchmark/target_benchmark.f90')]+[str(baseline/'lib'/lib) for lib in ['libaccelnet_target.a','libaccelnet.a','libAccelNetDescriptors.a']]+['-o',str(baseline/'accelnet-target-benchmark')]
 # The descriptor archive uses a lower-case output name on these builds.

 with (out/f'baseline-driver-{name}.log').open('w') as f:p=subprocess.run(cmd,cwd=mods,stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:raise SystemExit('FAIL baseline driver '+name)
