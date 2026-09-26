from pathlib import Path
import shutil,json,hashlib,difflib
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone');out=Path('/tmp/accelnet-angular-contraction');archive=root/'docs/validation/behler-contraction-2026-09-26';a=archive/'candidate-a';a.mkdir(exist_ok=True)
for p in out.glob('final-*'):
 if p.is_file():shutil.copy2(p,a/p.name)
for label,exe in [('gnu',Path('/tmp/accelnet-other-descriptors/build-gnu-serial/bin/accelnet-target-benchmark')),('nvhpc',Path('/tmp/accelnet-other-descriptors/build-nvhpc-serial/bin/accelnet-target-benchmark')),('gpu',Path('/tmp/accelnet-gpu-research/build-nvhpc/bin/accelnet-target-benchmark'))]:
 dst=out/('candidate-a-'+label);dst.mkdir(exist_ok=True);shutil.copy2(exe,dst/exe.name)
shutil.copytree(out/'source-stable/AccelNetPredictor/src',a/'src',dirs_exist_ok=True)
for name in ['accelnet_batch_target.f90','accelnet_target_descriptors.f90','accelnet_target_kernels.f90','accelnet_target_math.f90']:
 shutil.copy2(out/'source-reuse/AccelNetPredictor/src'/name,root/'AccelNetPredictor/src'/name)
shutil.copy2(out/'source-reuse/AccelNetPredictor/test/test_batch_target.f90',root/'AccelNetPredictor/test/test_batch_target.f90')
for name in ['reuse_gradients.py','check_reuse.py','reuse-progress.log','repeat-chebyshev.log']:
 shutil.copy2(out/name,archive/name)
shutil.copytree(out/'chebyshev-repeat-gnu',archive/'chebyshev-repeat-gnu',dirs_exist_ok=True)
# Reuse the fully defined baseline once; these preserved libraries remain revision 1.1.
s=(out/'build_final.py').read_text();s=s[:s.index(' baseline=')]
(out/'build_adopted.py').write_text(s)
s=(out/'final_checks.py').read_text()
# The full prior Chebyshev spot matrix plus the two-round follow-up are archived.
# Repeat the current implementation's checks once; keep adopted logs distinct.
s=s.replace("out=Path('/tmp/accelnet-angular-contraction')", "out=Path('/tmp/accelnet-angular-contraction/adopted-checks');out.mkdir(exist_ok=True)")
s=s.replace("before=out/f'baseline-{compiler}/accelnet-target-benchmark'", "before=Path('/tmp/accelnet-angular-contraction')/f'baseline-{compiler}/accelnet-target-benchmark'")
(out/'adopted_checks.py').write_text(s)
s=(out/'compare_final.py').read_text()
s=s.replace("for variant,exe in [('before',before),('after',after)]:", "for variant,exe in [('after',after)]:")
s=s.replace("f'final-{label}","f'adopted-{label}").replace("'final-timings.json'","'adopted-timings.json'")
(out/'compare_adopted.py').write_text(s)
