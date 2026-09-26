from pathlib import Path
import shutil
base=Path('/tmp/accelnet-angular-contraction');dst=base/'source-final';shutil.copytree(base/'source-pairs-local',dst,dirs_exist_ok=True)
root=dst/'AccelNetPredictor/src'
p=root/'accelnet_target_kernels.f90';s=p.read_text();a=s.index('    subroutine generic_forces');pre=s[:a];s=s[a:]
s=s.replace('        !$omp target teams distribute parallel do device(device)', '''        logical :: pair_once
        ! CPU serial builds benefit from evaluating both sides once. GPU edge
        ! ownership avoids contended force atomics; all scalar formulas are shared.
        pair_once = device == omp_get_initial_device()
        !$omp target teams distribute parallel do device(device)''',1)
s=s.replace('        !$omp& map(alloc:meta,nodes', '        !$omp& firstprivate(pair_once) map(alloc:meta,nodes')
s=s.replace('do k = j+1,offsets(row+1)-1', 'do k = merge(j+1,offsets(row),pair_once),offsets(row+1)-1\n                        if (k == j) cycle')
s=s.replace('hcoeff,degree,gradient,gradient_k)', 'hcoeff,degree,pair_once,gradient,gradient_k)')
s=s.replace('''                        do c=1,3
                            !$omp atomic update
                            edge_force(c,k) = edge_force(c,k)+gradient_k(c)
                        end do''','''                        if (pair_once) then
                            do c=1,3
                                !$omp atomic update
                                edge_force(c,k) = edge_force(c,k)+gradient_k(c)
                            end do
                        end if''')
s=s.replace('''            do c=1,3
                !$omp atomic update
                edge_force(c,j) = edge_force(c,j)+f(c)
            end do''','''            if (pair_once) then
                do c=1,3
                    !$omp atomic update
                    edge_force(c,j) = edge_force(c,j)+f(c)
                end do
            else
                edge_force(:,j) = f
            end if''');p.write_text(pre+s)
p=root/'accelnet_target_math.f90';s=p.read_text();a=s.index('    pure subroutine generic_pair_contracted');pre=s[:a];s=s[a:]
s=s.replace('coeff,degree,gradient,gradient_k)', 'coeff,degree,both,gradient,gradient_k)')
s=s.replace('        integer, intent(in) :: f(:),degree','        logical, intent(in) :: both\n        integer, intent(in) :: f(:),degree')
s=s.replace('        gradient_k=2*','        if (both) gradient_k=2*');p.write_text(pre+s)
# Clean up unused older pair helpers; their old source is archived.
p=root/'accelnet_target_math.f90';s=p.read_text()
import re
for name in ['generic_pair','generic_pair_cached','generic_pair_value']:
 s=re.sub(r'    pure (?:subroutine|real\(real64\) function) '+name+r'\([^\n]*(?:\n.*?)*?    end (?:subroutine|function)\n', '',s)
s=s.replace('public :: generic_radial, generic_radial_value, generic_pair, generic_pair_value, &\n            generic_pair_cached, generic_lj','public :: generic_radial, generic_radial_value, generic_lj')
p.write_text(s)
p=root/'accelnet_target_descriptors.f90';s=p.read_text().replace('! Real fields:', '! Fields 9..12: angular group, representative, next member, maximum polynomial degree.\n! Real fields:');p.write_text(s)
p=root/'accelnet_target_kernels.f90';s=p.read_text().replace('    subroutine generic_values', '''    ! Share pair geometry/radials within each angular group. Integer power sums
    ! remain local to the team, then each descriptor is written once. Stripping
    ! OpenMP directives gives the identical serial numerical loop.
    subroutine generic_values''',1).replace('    subroutine generic_contract', '''    ! Backpropagation coefficients are aggregated before visiting any force pair.
    ! The coefficient storage is O(centers * groups * degree), independent of edges.
    subroutine generic_contract''',1);p.write_text(s)
for p in root.glob('*.f90'):
 lines=[]
 for line in p.read_text().splitlines():
  while len(line)>130 and (not line.lstrip().startswith('!') or line.lstrip().startswith('!$omp')):
   k=line.rfind(',',0,110)
   if k<0:break
   lines.append(line[:k+1]+' &');line=('        !$omp& ' if line.lstrip().startswith('!$omp') else '            ')+line[k+1:].lstrip()
  lines.append(line)
 p.write_text('\n'.join(lines)+'\n')
