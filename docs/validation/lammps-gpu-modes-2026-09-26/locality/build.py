from pathlib import Path
import shlex,subprocess,shutil,hashlib,json
root=Path('/home/nagai/AccelNetGPU/AccelNet-clone')
base=Path('/tmp/accelnet-moment-analysis')
build=Path('/tmp/accelnet-gpu-research/build-nvhpc')
nv=Path('/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers/bin')
source=(root/'AccelNetPredictor/src/accelnet_target_kernels.f90').read_text()
start=source.index('            !$omp target teams distribute parallel do collapse(2)')
end=source.index('            !$omp end target teams distribute parallel do',start)
old=source[start:end]
link=shlex.split(Path('/tmp/accelnet-lammps-gpu/build/CMakeFiles/lmp.dir/link.txt').read_text())
records=[]
for tile in [1,32,128,512]:
 d=base/f'tile{tile}';d.mkdir(exist_ok=True)
 code=source.replace('integer :: entry, q, ax, ay, az, angular_neighbors, nm','integer :: entry, q, ax, ay, az, angular_neighbors, nm, tile, lane')
 new=old.replace('collapse(2)','collapse(3)').replace('private(s,j,ax,ay,az,mono,fcj,sj,self0,self1)','private(row,s,j,ax,ay,az,mono,fcj,sj,self0,self1)')
 new=new.replace('            do entry = 1, size(mp,2)\n                do row = 1, nrw',f'            do tile = 0, (nrw+{tile}-1)/{tile}-1\n              do entry = 1, size(mp,2)\n                do lane = 1, {tile}\n                    row = tile*{tile}+lane\n                    if (row > nrw) cycle')
 new += '            end do\n'
 code=code.replace(old,new)
 (d/'accelnet_target_kernels.f90').write_text(code)
 includes=[build/'AccelNetPredictor/target-modules',build/'AccelNetPredictor/modules',build/'AccelNetDescriptors/modules']
 command=[str(nv/'nvfortran'),'-fast','-O3','-mp=gpu','-gpu=cc90,cc120','-module',str(d),*[f'-I{x}' for x in includes],'-c',str(d/'accelnet_target_kernels.f90'),'-o',str(d/'accelnet_target_kernels.f90.o')]
 with (d/'compile.log').open('w') as f:subprocess.run(command,cwd=d,stdout=f,stderr=subprocess.STDOUT,check=True)
 lib=d/'libaccelnet_target.a';shutil.copy2(build/'lib/libaccelnet_target.a',lib)
 subprocess.run(['ar','r',str(lib),str(d/'accelnet_target_kernels.f90.o')],check=True)
 cmd=[str(lib) if x==str(build/'lib/libaccelnet_target.a') else x for x in link]
 cmd[cmd.index('-o')+1]=str(d/'lmp')
 with (d/'link.log').open('w') as f:subprocess.run(cmd,cwd='/tmp/accelnet-lammps-gpu/build',stdout=f,stderr=subprocess.STDOUT,check=True)
 records.append(dict(tile=tile,source_sha256=hashlib.sha256(code.encode()).hexdigest(),executable_sha256=hashlib.sha256((d/'lmp').read_bytes()).hexdigest()))
 print('built',tile,flush=True)
(base/'builds.json').write_text(json.dumps(records,indent=2)+'\n')
