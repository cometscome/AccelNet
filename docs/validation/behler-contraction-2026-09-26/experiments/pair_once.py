from pathlib import Path
import shutil
base=Path('/tmp/accelnet-angular-contraction'); shutil.copytree(base/'source',base/'source-pairs',dirs_exist_ok=True)
root=base/'source-pairs/AccelNetPredictor/src'
p=root/'accelnet_target_kernels.f90';s=p.read_text();s=s.replace('                    edge_row(j) = row','                    edge_row(j) = row\n                    edge_force(:,j) = 0',1)
start=s.index('    subroutine generic_forces');a=s[:start];f=s[start:]
f=f.replace('group,ag,degree','group,ag,degree,c').replace('gradient(3),v12','gradient(3),gradient_k(3),v12').replace('f,v,dv,gradient,v12','f,v,dv,gradient,gradient_k,v12')
f=f.replace('do k = offsets(row),offsets(row+1)-1','do k = j+1,offsets(row+1)-1').replace('                        if (k == j) cycle\n','')
f=f.replace('radial_cache(j,group,2),angular_coeff,row,ag,degree,gradient)', 'radial_cache(j,group,2),radial_cache(k,group,2),angular_coeff,row,ag,degree,gradient,gradient_k)')
f=f.replace('                        f = f+gradient','''                        f = f+gradient
                        do c=1,3
                            !$omp atomic update
                            edge_force(c,k) = edge_force(c,k)+gradient_k(c)
                        end do''')
f=f.replace('            edge_force(:,j) = f+radial_force*geom(1:3,j)','''            f = f+radial_force*geom(1:3,j)
            do c=1,3
                !$omp atomic update
                edge_force(c,j) = edge_force(c,j)+f(c)
            end do''')
p.write_text(a+f)
p=root/'accelnet_target_math.f90';s=p.read_text();a=s.index('    pure subroutine generic_pair_contracted');prefix=s[:a];s=s[a:]
s=s.replace('qj,qk,dqj,coeff,row,ag,degree,gradient)', 'qj,qk,dqj,dqk,coeff,row,ag,degree,gradient,gradient_k)')
s=s.replace('qj,qk,dqj,coeff(:,:,:)', 'qj,qk,dqj,dqk,coeff(:,:,:)').replace(':: gradient(3)',':: gradient(3),gradient_k(3)').replace('        gradient=0','        gradient=0; gradient_k=0')
s=s.replace('    end subroutine', '''        gradient_k=2*(da*product*(uj-cosine*uk)/rk+a*qj*(dqk*qjk*uk+qk*dqjk*ujk))
    end subroutine''')
p.write_text(prefix+s)
for p in root.glob('*.f90'):
 out=[]
 for line in p.read_text().splitlines():
  while len(line)>130 and (not line.lstrip().startswith('!') or line.lstrip().startswith('!$omp')):
   k=line.rfind(',',0,110)
   if k<0:break
   out.append(line[:k+1]+' &');line=('        !$omp& ' if line.lstrip().startswith('!$omp') else '            ')+line[k+1:].lstrip()
  out.append(line)
 p.write_text('\n'.join(out)+'\n')
