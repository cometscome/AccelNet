#!/usr/bin/env python3
"""Build a controlled CPU-only performance experiment without editing the library.

cached: cache the identical angular radial values/derivatives before pair loops.
lj-powers: additionally replace LJ integer powers with explicit multiplications.
lj-fused: instead evaluate the two LJ features together, sharing cutoff/powers.

Use only a build configured with ACCELNET_TARGET_SERIAL=ON. The output includes
replacement sources, a static library, the unchanged timing driver, and the
unchanged numerical/finite-difference test. No production files are changed.
"""
from pathlib import Path
import argparse
import re
import shlex
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--build', required=True, type=Path)
p.add_argument('--output', required=True, type=Path)
p.add_argument('--variant', choices=['cached', 'lj-powers', 'lj-fused'], default='cached')
a = p.parse_args()
root = Path(__file__).resolve().parents[2]
build = a.build.resolve()
out = a.output.resolve()
out.mkdir(parents=True, exist_ok=True)
cache = (build/'CMakeCache.txt').read_text()
if 'ACCELNET_TARGET_SERIAL:BOOL=ON' not in cache:
    p.error('This experiment requires ACCELNET_TARGET_SERIAL=ON (one serial CPU thread)')
match = re.search(r'^CMAKE_Fortran_COMPILER:[^=]+=(.+)$', cache, re.M)
if not match:
    p.error('Missing Fortran compiler in CMakeCache.txt')
compiler_path = match[1]
compiler = 'nvhpc' if 'nvfortran' in compiler_path else 'gnu'
lj_powers = a.variant == 'lj-powers'
fused_lj = a.variant == 'lj-fused'
# This historical source-rewriting experiment predates persistent group caches.
# Refuse a newer kernel ABI instead of producing a misleading comparison.
if 'radial_cache' in (root/'AccelNetPredictor/src/accelnet_target_kernels.f90').read_text():
    p.error('The common kernel already has persistent radial caching and LJ fusion. '
            'Use compare_target_backends.py for current code; historical diagnostic '
            'sources/results are in docs/validation/gpu-descriptors-cpu-comparison-2026-09-26/.')
math=(root/'AccelNetPredictor/src/accelnet_target_math.f90').read_text()
start=math.index('    pure subroutine generic_pair(');end=math.index('    end subroutine',start)+len('    end subroutine')
fn=math[start:end].replace('generic_pair(', 'generic_pair_cached(').replace('rj, rk, value, derivative)', 'rj, rk, qj, qk, dqj, value, derivative)')
fn=fn.replace('        real(real64), intent(out) :: value,derivative(3)','        real(real64), intent(in) :: qj,qk,dqj\n        real(real64), intent(out) :: value,derivative(3)')
fn=fn.replace(':: qj,qk,dqj,dqk,qjk', ':: qjk')
fn=fn.replace('        call generic_radial(2,rj,f(4),p,qj,dqj)\n','').replace('        call generic_radial(2,rk,f(4),p,qk,dqk)\n','')
math=math.replace('    public :: generic_radial, generic_pair','    public :: generic_radial, generic_pair, generic_pair_cached').replace('end module',fn+'\nend module')
kernel=(root/'AccelNetPredictor/src/accelnet_target_kernels.f90').read_text()
for name in ['generic_values','generic_forces']:
 start=kernel.index('    subroutine '+name);end=kernel.index('    end subroutine',start)
 sub=kernel[start:end]
 pos=sub.index('        !$omp target teams')
 declarations='        real(real64), allocatable :: qcache(:,:),dqcache(:,:)\n'
 if name=='generic_values':
  pre='''        allocate(qcache(offsets(nrw+1)-1,size(features,2)),dqcache(offsets(nrw+1)-1,size(features,2)))
       do row = 1,nrw
           s = species(centers(row))
           do b = 1,nodes(1,s)
               if (features(1,b,s) /= 4 .and. features(1,b,s) /= 5) cycle
               do j = offsets(row),offsets(row+1)-1
                   call generic_radial(2,geom(4,j),features(4,b,s),fp(:,b,s),qcache(j,b),dqcache(j,b))
               end do
           end do
       end do
'''
 else:
  pre='''        allocate(qcache(nedges,size(features,2)),dqcache(nedges,size(features,2)))
       do j = 1,nedges
           row = edge_row(j); s = species(centers(row))
           do b = 1,nodes(1,s)
               if (features(1,b,s) /= 4 .and. features(1,b,s) /= 5) cycle
               call generic_radial(2,geom(4,j),features(4,b,s),fp(:,b,s),qcache(j,b),dqcache(j,b))
           end do
       end do
'''
 sub=sub[:pos]+declarations+pre+sub[pos:]
 sub=sub.replace('call generic_pair(', 'call generic_pair_cached(').replace('geom(4,j),geom(4,k),v,gradient)', 'geom(4,j),geom(4,k),qcache(j,b),qcache(k,b),dqcache(j,b),v,gradient)')
 kernel=kernel[:start]+sub+kernel[end:]
if lj_powers:
 math=math.replace('h = r**(-6);', 'h = 1.0_real64/(r*r); h = h*h*h;').replace('h = r**(-12);', 'h = 1.0_real64/(r*r); h = h*h*h; h = h*h;')
if fused_lj:
 helper="""    pure subroutine generic_lj(r,ct,p,force,v,v12,dv,dv12)
       integer, intent(in) :: ct
       real(real64), intent(in) :: r,p(7)
       logical, intent(in) :: force
       real(real64), intent(out) :: v,v12,dv,dv12
       real(real64) :: ir2,ir6,ir12,fc,dfc
       v=0;v12=0;dv=0;dv12=0
       if (r <= 1e-12_real64 .or. r > p(1)) return
       ir2=1.0_real64/(r*r);ir6=ir2*ir2*ir2;ir12=ir6*ir6
       fc=target_cutoff_value(r,p(1),ct,p(7))
       v=fc*ir6;v12=fc*ir12
       if (force) then
           dfc=target_cutoff_derivative(r,p(1),ct,p(7))
           dv=ir6*(dfc-6*fc/r);dv12=ir12*(dfc-12*fc/r)
       end if
   end subroutine
"""
 math=math.replace('public :: generic_radial,', 'public :: generic_lj, generic_radial,').replace('end module',helper+'end module')
 start=kernel.index('    subroutine generic_values');end=kernel.index('    subroutine generic_forces',start)
 value=kernel[start:end]
 value=value.replace(':: total,v,dv,gradient(3)', ':: total,v,dv,gradient(3),total12,v12,dv12')
 value=value.replace('                total = 0','                if (kind == 7) cycle\n                total = 0; total12 = 0')
 value=value.replace('                        call generic_radial(kind,geom(4,j),features(4,b,s),fp(:,b,s),v,dv)',
     '                        if (kind == 6) then\n                            call generic_lj(geom(4,j),features(4,b,s),fp(:,b,s),.false.,v,v12,dv,dv12)\n                            total12=total12+v12\n                        else\n                            call generic_radial(kind,geom(4,j),features(4,b,s),fp(:,b,s),v,dv)\n                        end if')
 value=value.replace('                g(row,b) = total','                g(row,b) = total\n                if (kind == 6) g(row,b+1) = total12')
 kernel=kernel[:start]+value+kernel[end:]
 start=kernel.index('    subroutine generic_forces');force=kernel[start:]
 force=force.replace(':: f(3),v,dv,gradient(3)', ':: f(3),v,dv,gradient(3),v12,dv12')
 force=force.replace('                if (kind /= 4', '                if (kind == 7) cycle\n                if (kind /= 4')
 force=force.replace('                    call generic_radial(kind,geom(4,j),features(4,b,s),fp(:,b,s),v,dv)\n                    f = f+g(row,b)*dv*geom(1:3,j)',
     '                    if (kind == 6) then\n                        call generic_lj(geom(4,j),features(4,b,s),fp(:,b,s),.true.,v,v12,dv,dv12)\n                        f=f+(g(row,b)*dv+g(row,b+1)*dv12)*geom(1:3,j)\n                    else\n                        call generic_radial(kind,geom(4,j),features(4,b,s),fp(:,b,s),v,dv)\n                        f=f+g(row,b)*dv*geom(1:3,j)\n                    end if')
 kernel=kernel[:start]+force
(out/'accelnet_target_math.f90').write_text(math)
(out/'accelnet_target_kernels.f90').write_text(kernel)
cc = compiler_path
flags=['-O3','-ffree-line-length-none','-J',str(out)] if compiler=='gnu' else ['-fast','-O3','-module',str(out)]
includes=[out,build/'AccelNetPredictor/target-modules',build/'AccelNetPredictor/modules',build/'AccelNetDescriptors/modules',root/'AccelNetDescriptors/src/shared']
for src in ['accelnet_target_math.f90','accelnet_target_kernels.f90']:
 cmd=[cc,*flags,*sum((['-I',str(d)] for d in includes),[]),'-c',str(out/src),'-o',str(out/(src+'.o'))]
 subprocess.run(cmd,check=True)
import shutil
lib=out/'libaccelnet_target_cached.a';shutil.copyfile(build/'lib/libaccelnet_target.a',lib)
subprocess.run(['ar','r',str(lib),str(out/'accelnet_target_math.f90.o'),str(out/'accelnet_target_kernels.f90.o')],check=True)
for target in ['accelnet-target-benchmark','test_batch_target']:
 args=shlex.split((build/f'AccelNetPredictor/CMakeFiles/{target}.dir/link.txt').read_text())
 args[args.index('-o')+1]=str(out/target)
 args=[str(lib) if x.endswith('/libaccelnet_target.a') else x for x in args]
 subprocess.run(args,cwd=build/'AccelNetPredictor',check=True)
print('Built cache experiment',compiler,flush=True)
