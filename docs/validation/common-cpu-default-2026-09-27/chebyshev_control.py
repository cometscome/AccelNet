from pathlib import Path
import subprocess,sys
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-common-default/chebyshev-control')
for compiler in ['gnu','nvhpc']:
 cmd=[sys.executable,str(root/'AccelNetPredictor/benchmark/compare_chebyshev_variants.py'),'--variant','before','cpu-shared',f'/tmp/accelnet-common-default/before-{compiler}','--variant','after','cpu-shared',f'/tmp/accelnet-other-descriptors/build-{compiler}-serial/bin/accelnet-target-benchmark','--family','chebyshev','--sizes','512','4096','--orders','8','--modes','1','--rounds','2','--seconds','0.10','--cpu','6','--output',str(out/compiler)]
 subprocess.run(cmd,cwd=root,check=True)
