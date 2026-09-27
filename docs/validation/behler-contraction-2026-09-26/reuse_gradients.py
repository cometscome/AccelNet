from pathlib import Path
import shutil,re
base=Path('/tmp/accelnet-angular-contraction');dst=base/'source-reuse';shutil.copytree(base/'source-stable',dst,dirs_exist_ok=True)
root=dst/'AccelNetPredictor/src';before=Path('/tmp/accelnet-jacobian-trial/before-common/src')
p=root/'accelnet_batch_target.f90';s=(before/p.name).read_text().replace('self%features(8,mn,ns)','self%features(14,mn,ns)');p.write_text(s)
p=root/'accelnet_target_descriptors.f90';s=p.read_text().replace('fi(12,n)','fi(14,n)')
s=s.replace('fi(5,b) > 0 .and. fi(5,b) <= 16', 'fi(5,b) > 0 .and. fi(5,b) <= 16 .and. fr(5,b) == real(fi(5,b),real64)')
s=s.replace("        status = 0; message = ''",'''        ! Equal-power coefficients occupy their original first NN-gradient slot.
        ! Fields 13/14 link only equal powers, while field 11 links the full basis.
        do b = 1,n
            if (fi(10,b) == 0) cycle
            k = fi(10,b)
            do while (k /= b)
                if (fr(5,k) == fr(5,b)) exit
                k = fi(11,k)
            end do
            fi(13,b) = k
            if (k == b) cycle
            do while (fi(14,k) /= 0)
                k = fi(14,k)
            end do
            fi(14,k) = b
        end do
        status = 0; message = '' ''')
s=s.replace('! Real fields:', '! Fields 13/14: equal-power representative and next equal-power member.\n! Real fields:');p.write_text(s)
p=root/'accelnet_target_kernels.f90';s=p.read_text();a=s.index('    subroutine generic_values');pre=s[:a];tail=s[a:]
pre=pre.replace('radial_cache, angular_coeff, use_moment','radial_cache, use_moment').replace(', angular_coeff(:,:,:)','')
pre=pre.replace('nm, ag, degree, bb, d','nm, bb').replace('self0,self1,ag,degree,bb,d)', 'self0,self1,bb)')
a=pre.index('            ! Fuse angular coefficient aggregation');end=pre.index('            nr = meta(1,s); na = meta(2,s); multi = meta(3,s)',a)
pre=pre[:a]+'''            ! Reuse NN gradients for equal-power coefficient sums. Distinct
            ! powers already have the desired coefficient, so leave them in place.
            ! No additional array, allocation, transfer, or kernel launch.
            if (meta(10,s) == 1) then
                do b=1,dim
                    if (features(13,b,s) /= b .or. features(14,b,s) == 0) cycle
                    v=g(row,b); bb=features(14,b,s)
                    do while (bb /= 0)
                        v=v+g(row,bb)
                        bb=features(14,bb,s)
                    end do
                    g(row,b)=v
                end do
            end if
'''+pre[end:]
pre=pre.replace('radial_cache,angular_coeff,g,edge_force)', 'radial_cache,g,edge_force)').replace('                    edge_force(:,j) = 0\n','')
tail=tail.replace('radial_cache,angular_coeff,g,edge_force)', 'radial_cache,g,edge_force)').replace(',angular_coeff(:,:,:)','')
tail=tail.replace('group,ag,degree,d,c','group,degree,d,c,bb')
tail=tail.replace('pair_once = device == omp_get_initial_device()', '''pair_once = device == omp_get_initial_device() .and. any(features(10,:,:) > 0)
        if (pair_once) then
            do j=1,nedges
                row=edge_row(j); s=species(centers(row))
                if (meta(10,s) == 1) edge_force(:,j)=0
            end do
        end if''')
tail=tail.replace('''                    ag = features(9,b,s); degree = features(12,b,s)
                    do d=0,degree
                        hcoeff(d) = angular_coeff(row,ag,d+1)
                    end do''','''                    degree = features(12,b,s)
                    hcoeff(:degree) = 0
                    if (degree > 0) then
                        bb=b
                        do while (bb /= 0)
                            if (features(13,bb,s) == bb) hcoeff(features(5,bb,s))=g(row,bb)
                            bb=features(11,bb,s)
                        end do
                    else
                        hcoeff(0)=g(row,b)
                    end if''')
s=pre+tail;assert 'angular_coeff' not in s;p.write_text(s)
p=dst/'AccelNetPredictor/test/test_batch_target.f90';s=p.read_text().replace("work%allocations() == nw+1,'angular coefficient degree growth'", "work%allocations() == nw,'angular coefficient storage reuse'");p.write_text(s)
for p in root.glob('*.f90'):
 lines=[]
 for line in p.read_text().splitlines():
  while len(line)>130 and (not line.lstrip().startswith('!') or line.lstrip().startswith('!$omp')):
   k=line.rfind(',',0,110)
   if k<0:break
   lines.append(line[:k+1]+' &');line=('        !$omp& ' if line.lstrip().startswith('!$omp') else '            ')+line[k+1:].lstrip()
  lines.append(line)
 p.write_text('\n'.join(lines)+'\n')
