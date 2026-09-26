from pathlib import Path
import shutil,re
base=Path('/tmp/accelnet-angular-contraction');dst=base/'source-stable';shutil.copytree(base/'source-final',dst,dirs_exist_ok=True)
root=dst/'AccelNetPredictor/src';before=Path('/tmp/accelnet-jacobian-trial/before-common/src')
p=root/'accelnet_target_descriptors.f90';s=p.read_text();s=s.replace("        status = 0; message = ''",'''        ! A single distinct power needs no polynomial recurrence. Keep the
        ! original power helper, and only aggregate its scalar NN coefficient.
        do b = 1,n
            if (fi(10,b) /= b .or. fi(12,b) == 0) cycle
            k = fi(11,b)
            do while (k /= 0)
                if (fi(5,k) /= fi(5,b)) exit
                k = fi(11,k)
            end do
            if (k == 0) fi(12,b) = 0
        end do
        status = 0; message = '' ''');p.write_text(s)
p=root/'accelnet_target_kernels.f90';s=p.read_text()
s=s.replace('integer :: entry, q, ax, ay, az, angular_neighbors, nm', 'integer :: entry, q, ax, ay, az, angular_neighbors, nm, ag, degree, bb, d')
s=s.replace('private(s,dim,l,nin,nout,o,j,i,z,v,b,q,nr,na,multi,entry,self0,self1)', 'private(s,dim,l,nin,nout,o,j,i,z,v,b,q,nr,na,multi,entry,self0,self1,ag,degree,bb,d)')
s=s.replace('''        if (any(meta(10,:) == 1)) call generic_values(device,nrw,meta,nodes,features,feature_params, &
            local_species,species,centers,offsets,indices,geom,radial_cache,g)''','''        if (any(meta(10,:) == 1)) then
            if (any(features(12,:,:) > 0)) then
                call generic_values_grouped(device,nrw,meta,nodes,features,feature_params, &
                    local_species,species,centers,offsets,indices,geom,radial_cache,g)
            else
                call generic_values(device,nrw,meta,nodes,features,feature_params, &
                    local_species,species,centers,offsets,indices,geom,radial_cache,g)
            end if
        end if''')
needle='''            nr = meta(1,s); na = meta(2,s); multi = meta(3,s)
            if (use_moment(row) == 1) then'''
replacement='''            ! Fuse angular coefficient aggregation into the existing NN launch,
            ! as for Chebyshev moments. No additional target launch or transfers.
            if (meta(10,s) == 1) then
                do b=1,dim
                    if (features(10,b,s) /= b) cycle
                    ag=features(9,b,s); degree=features(12,b,s)
                    do d=1,degree+1
                        angular_coeff(row,ag,d)=0
                    end do
                    bb=b
                    do while (bb /= 0)
                        d=0
                        if (degree > 0) d=features(5,bb,s)
                        angular_coeff(row,ag,d+1)=angular_coeff(row,ag,d+1)+g(row,bb)
                        bb=features(11,bb,s)
                    end do
                end do
            end if
'''+needle
assert needle in s;s=s.replace(needle,replacement)
s=s.replace('        if (any(meta(10,:) == 1)) call generic_contract(device,nrw,features,species,centers,g,angular_coeff)\n','')
a=s.index('    ! Backpropagation coefficients are aggregated');e=s.index('    subroutine generic_forces',a);s=s[:a]+s[e:]
s=s.replace('    subroutine generic_values(', '    subroutine generic_values_grouped(',1)
old=(before/'accelnet_target_kernels.f90').read_text();a=old.index('    subroutine generic_values');e=old.index('    subroutine generic_forces',a)
s=s.replace('    subroutine generic_values_grouped',old[a:e]+'    subroutine generic_values_grouped',1);p.write_text(s)
p=root/'accelnet_target_math.f90';s=p.read_text();old=(before/p.name).read_text();a=old.index('    pure real(real64) function generic_pair_value');e=old.index('    end function',a)+len('    end function\n')
s=s.replace('public :: generic_pair_geometry, generic_pair_contracted','public :: generic_pair_geometry, generic_pair_contracted, generic_pair_value')
s=s.replace('end module',old[a:e]+'end module');p.write_text(s)
# Dense, non-duplicate basis benchmark uses ordinary model APIs and the same NN builder.
p=dst/'AccelNetPredictor/test/batch_test_support.f90';s=p.read_text().replace('degree, v','degree, v, pair, sign, t1, t2')
needle="            case ('behler', 'g4', 'g5', 'lj-behler')"
new='''            case ('g4-series', 'g5-series')
                call initialize_behler_config(behler,2)
                do pair = 1,3
                    t1 = merge(2,1,pair == 3); t2 = merge(1,2,pair == 1)
                    do sign = 1,2
                        do j = 1,max(2,degree)
                            if (family == 'g4-series') then
                                call add_g4(behler,t1,t2,3.4_real64,real(2*sign-3,real64), &
                                    real(j,real64),0.2_real64,0.4_real64)
                            else
                                call add_g5(behler,t1,t2,3.4_real64,real(2*sign-3,real64), &
                                    real(j,real64),0.2_real64,0.4_real64)
                            end if
                        end do
                    end do
                end do
                call add_behler(model%setups(s)%model,behler)
'''+needle
assert needle in s;s=s.replace(needle,new);p.write_text(s)
for p in root.glob('*.f90'):
 lines=[]
 for line in p.read_text().splitlines():
  while len(line)>130 and (not line.lstrip().startswith('!') or line.lstrip().startswith('!$omp')):
   k=line.rfind(',',0,110)
   if k<0:break
   lines.append(line[:k+1]+' &');line=('        !$omp& ' if line.lstrip().startswith('!$omp') else '            ')+line[k+1:].lstrip()
  lines.append(line)
 p.write_text('\n'.join(lines)+'\n')
