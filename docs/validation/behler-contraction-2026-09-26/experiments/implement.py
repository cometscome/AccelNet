from pathlib import Path
root=Path('/tmp/accelnet-angular-contraction/source/AccelNetPredictor/src')
p=root/'accelnet_target_descriptors.f90';s=p.read_text().replace('groups,k','groups,k,r,last').replace('fi(8,n)','fi(12,n)')
s=s.replace("        status = 0; message = ''",'''        ! Contract only features with identical non-angular factors and lambda.
        ! Positive integer zeta <= 16 shares a polynomial in t=(1+lambda*c)/2.
        ! Fractional/high zeta retains its exact power and only identical powers combine.
        groups = 0
        do b = 1,n
            if (fi(7,b) == 0) cycle
            r = 0
            do k = 1,b-1
                if (fi(10,k) /= k) cycle
                if (fi(1,b) /= fi(1,k) .or. fi(7,b) /= fi(7,k)) cycle
                if (minval(fi(2:3,b)) /= minval(fi(2:3,k))) cycle
                if (maxval(fi(2:3,b)) /= maxval(fi(2:3,k))) cycle
                if (fr(4,b) /= fr(4,k)) cycle
                if (fi(5,b) > 0 .and. fi(5,b) <= 16 .and. fi(12,k) > 0) then
                    r = k
                else if (fi(12,k) == 0 .and. fr(5,b) == fr(5,k)) then
                    r = k
                end if
                if (r /= 0) exit
            end do
            if (r == 0) then
                groups = groups+1
                fi(9,b) = groups; fi(10,b) = b
                if (fi(5,b) > 0 .and. fi(5,b) <= 16) fi(12,b) = fi(5,b)
            else
                fi(9:10,b) = fi(9:10,r)
                last = r
                do while (fi(11,last) /= 0)
                    last = fi(11,last)
                end do
                fi(11,last) = b
                if (fi(12,r) > 0) fi(12,r) = max(fi(12,r),fi(5,b))
            end if
        end do
        status = 0; message = '' ''')
p.write_text(s)
p=root/'accelnet_batch_target.f90';s=p.read_text().replace('self%features(8,mn,ns)','self%features(12,mn,ns)')
s=s.replace('powers(:,:,:), radial(:,:,:)','powers(:,:,:), radial(:,:,:), angular_coeff(:,:,:)')
s=s.replace('self%powers,self%radial,','self%powers,self%radial,self%angular_coeff,').replace('work%powers,work%radial,','work%powers,work%radial,work%angular_coeff,').replace('powers,radial,species','powers,radial,angular_coeff,species')
s=s.replace('integer :: n,mn,ml,ne,na,mm,mpower,ngroups','integer :: n,mn,ml,ne,na,mm,mpower,ngroups,nag,nad')
s=s.replace('        grow = .true.', '''        nag = max(1,maxval(model%features(9,:,:)))
        nad = 1+maxval(model%features(12,:,:))
        grow = .true.''')
s=s.replace('            ngroups = max(ngroups,size(work%radial,2))', '''            ngroups = max(ngroups,size(work%radial,2))
            nag = max(nag,size(work%angular_coeff,2)); nad = max(nad,size(work%angular_coeff,3))''')
s=s.replace(' .or. ngroups > size(work%radial,2)', ' .or. ngroups > size(work%radial,2) .or. &\n                nag > size(work%angular_coeff,2) .or. nad > size(work%angular_coeff,3)')
s=s.replace('work%radial(ne,ngroups,2),','work%radial(ne,ngroups,2),work%angular_coeff(n,nag,nad),')
p.write_text(s)
p=root/'accelnet_target_kernels.f90';s=p.read_text()
s=s.replace('powers, radial_cache, use_moment','powers, radial_cache, angular_coeff, use_moment')
s=s.replace('powers(:,:,:), radial_cache(:,:,:)','powers(:,:,:), radial_cache(:,:,:), angular_coeff(:,:,:)')
s=s.replace('        if (any(meta(10,:) == 1)) call generic_forces', '''        if (any(meta(10,:) == 1)) call generic_contract(device,nrw,features,species,centers,g,angular_coeff)
        if (any(meta(10,:) == 1)) call generic_forces''')
s=s.replace('geom,edge_row,radial_cache,g,edge_force)', 'geom,edge_row,radial_cache,angular_coeff,g,edge_force)')
start=s.index('    subroutine generic_values');stop=s.index('    subroutine generic_forces',start)
values=s[start:stop]; force=s[stop:]
values=values.replace('kind,t1,t2,group','kind,t1,t2,group,bb').replace('gradient(3),total12','gradient(3),qjk,dqjk,ujk(3),cosine,total12')
values=values.replace('gradient,total12','gradient,qjk,dqjk,ujk,cosine,total12')
values=values.replace('''                if (kind == 7) cycle''','''                if (kind == 7) cycle
                if (features(10,b,s) /= 0 .and. features(10,b,s) /= b) cycle
                if (features(10,b,s) == b) then
                    bb = b
                    do while (bb /= 0)
                        g(row,bb) = 0
                        bb = features(11,bb,s)
                    end do
                end if''')
a=values.index('                            v = generic_pair_value');z=values.index('                            total = total+v',a)+len('                            total = total+v')
values=values[:a]+'''                            call generic_pair_geometry(features(:,b,s),fp(:,b,s),geom(1:3,j),geom(1:3,k), &
                                geom(4,j),geom(4,k),cosine,qjk,dqjk,ujk,.false.)
                            if (qjk == 0) cycle
                            v = 2*radial_cache(j,group,1)*radial_cache(k,group,1)*qjk
                            bb = b
                            do while (bb /= 0)
                                g(row,bb) = g(row,bb)+v*angular_value(cosine,fp(4,bb,s),fp(5,bb,s),features(5,bb,s))
                                bb = features(11,bb,s)
                            end do'''+values[z:]
values=values.replace('                g(row,b) = total','                if (features(10,b,s) == 0) g(row,b) = total')
force=force.replace('fp(:,:,:),radial_cache(:,:,:),geom(:,:),g(:,:)', 'fp(:,:,:),radial_cache(:,:,:),angular_coeff(:,:,:),geom(:,:),g(:,:)')
force=force.replace('kind,t1,t2,group','kind,t1,t2,group,ag,degree')
force=force.replace('                if (kind == 7) cycle', '''                if (kind == 7) cycle
                if (features(10,b,s) /= 0 .and. features(10,b,s) /= b) cycle''')
force=force.replace('''                    if (tj /= t1 .and. tj /= t2) cycle''','''                    if (tj /= t1 .and. tj /= t2) cycle
                    ag = features(9,b,s); degree = features(12,b,s)''')
force=force.replace('call generic_pair_cached(', 'call generic_pair_contracted(')
force=force.replace('radial_cache(j,group,2),gradient)', 'radial_cache(j,group,2),angular_coeff,row,ag,degree,gradient)')
force=force.replace('f = f+g(row,b)*gradient','f = f+gradient')
contract='''    subroutine generic_contract(device,nrw,features,species,centers,g,coeff)
        integer, intent(in) :: device,nrw
        integer, contiguous, intent(in) :: features(:,:,:),species(:),centers(:)
        real(real64), contiguous, intent(in) :: g(:,:)
        real(real64), contiguous, intent(inout) :: coeff(:,:,:)
        integer :: row,b,bb,s,ag,d,n
        !$omp target teams distribute parallel do collapse(2) device(device) if(device /= omp_get_initial_device()) &
        !$omp& map(alloc:features,species,centers,g,coeff) private(bb,s,ag,d,n)
        do row=1,nrw
            do b=1,size(features,2)
                s=species(centers(row))
                if (features(10,b,s) /= b) cycle
                ag=features(9,b,s); n=features(12,b,s)
                do d=1,n+1
                    coeff(row,ag,d)=0
                end do
                bb=b
                do while (bb /= 0)
                    d=0
                    if (n > 0) d=features(5,bb,s)
                    coeff(row,ag,d+1)=coeff(row,ag,d+1)+g(row,bb)
                    bb=features(11,bb,s)
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine

'''
s=s[:start]+values+contract+force;p.write_text(s)
p=root/'accelnet_target_math.f90';s=p.read_text().replace('public :: generic_radial,','public :: generic_pair_geometry, generic_pair_contracted\n    public :: generic_radial,')
helpers='''
    pure subroutine generic_pair_geometry(f,p,uj,uk,rj,rk,cosine,qjk,dqjk,ujk,force)
        !$omp declare target
        integer, intent(in) :: f(:)
        real(real64), intent(in) :: p(7),uj(3),uk(3),rj,rk
        logical, intent(in) :: force
        real(real64), intent(out) :: cosine,qjk,dqjk,ujk(3)
        real(real64) :: rjk
        cosine=0; qjk=0; dqjk=0; ujk=0
        if (rj <= 1e-12_real64 .or. rk <= 1e-12_real64 .or. rj > p(1) .or. rk > p(1)) return
        qjk=1
        if (f(1) == 4) then
            ujk=rk*uk-rj*uj; rjk=sqrt(sum(ujk**2))
            if (rjk <= 1e-12_real64 .or. rjk > p(1)) then
                qjk=0; return
            end if
            if (force) then
                ujk=ujk/rjk
                call generic_radial(2,rjk,f(4),p,qjk,dqjk)
            else
                qjk=generic_radial_value(2,rjk,f(4),p)
            end if
        end if
        cosine=max(-1.0_real64,min(1.0_real64,sum(uj*uk)))
    end subroutine

    pure subroutine generic_pair_contracted(f,p,uj,uk,rj,rk,qj,qk,dqj,coeff,row,ag,degree,gradient)
        !$omp declare target
        integer, intent(in) :: f(:),row,ag,degree
        real(real64), intent(in) :: p(7),uj(3),uk(3),rj,rk,qj,qk,dqj,coeff(:,:,:)
        real(real64), intent(out) :: gradient(3)
        real(real64) :: cosine,qjk,dqjk,ujk(3),a,da,t,product
        integer :: d
        gradient=0
        call generic_pair_geometry(f,p,uj,uk,rj,rk,cosine,qjk,dqjk,ujk,.true.)
        if (qjk == 0) return
        if (degree > 0) then
            t=0.5_real64*(1+p(4)*cosine)
            a=coeff(row,ag,degree+1); da=0
            do d=degree-1,0,-1
                da=da*t+a
                a=a*t+coeff(row,ag,d+1)
            end do
            da=da*0.5_real64*p(4)
        else
            call angular_power(cosine,p(4),p(5),f(5),0.5_real64*p(5)*p(4),a,da)
            a=a*coeff(row,ag,1); da=da*coeff(row,ag,1)
        end if
        product=qj*qk*qjk
        gradient=2*(da*product*(uk-cosine*uj)/rj+a*qk*(dqj*qjk*uj-qj*dqjk*ujk))
    end subroutine
'''
s=s.replace('end module',helpers+'end module');p.write_text(s)
for p in root.glob('*.f90'):
 out=[]
 for line in p.read_text().splitlines():
  while len(line)>130 and (not line.lstrip().startswith('!') or line.lstrip().startswith('!$omp')):
   k=line.rfind(',',0,110)
   if k<0:break
   out.append(line[:k+1]+' &');line=('        !$omp& ' if line.lstrip().startswith('!$omp') else '            ')+line[k+1:].lstrip()
  out.append(line)
 p.write_text('\n'.join(out)+'\n')
