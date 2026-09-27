from pathlib import Path
root=Path('/tmp/accelnet-jacobian-trial/source/AccelNetPredictor/src')
p=root/'accelnet_batch_target.f90';s=p.read_text()
s=s.replace('powers(:,:,:), radial(:,:,:)','powers(:,:,:), radial(:,:,:), jacobian(:,:,:)')
s=s.replace('self%powers,self%radial,','self%powers,self%radial,self%jacobian,').replace('work%powers,work%radial,','work%powers,work%radial,work%jacobian,')
s=s.replace('powers,radial,species','powers,radial,jacobian,species')
s=s.replace('integer :: n,mn,ml,ne,na,mm,mpower,ngroups','integer :: n,mn,ml,ne,na,mm,mpower,ngroups,njf,nje')
s=s.replace('        grow = .true.','''        njf = max(1,maxval(model%nodes(1,:),mask=model%meta(10,:) == 1))
        nje = 1
        if (any(model%meta(10,:) == 1)) nje = ne
        grow = .true.''')
s=s.replace('            ngroups = max(ngroups,size(work%radial,2))','''            ngroups = max(ngroups,size(work%radial,2))
            njf = max(njf,size(work%jacobian,2)); nje = max(nje,size(work%jacobian,1))''')
s=s.replace(' .or. ngroups > size(work%radial,2)',' .or. ngroups > size(work%radial,2) .or. &\n                njf > size(work%jacobian,2) .or. nje > size(work%jacobian,1)')
s=s.replace('work%radial(ne,ngroups,2),','work%radial(ne,ngroups,2),work%jacobian(nje,njf,3),')
p.write_text(s)
p=root/'accelnet_target_math.f90';s=p.read_text().replace('public :: generic_radial,','public :: generic_pair_both, generic_radial,')
helper='''
    pure subroutine generic_pair_both(f,p,uj,uk,rj,rk,qj,qk,dqj,dqk,value,dj,dk)
        !$omp declare target
        integer, intent(in) :: f(:)
        real(real64), intent(in) :: p(7),uj(3),uk(3),rj,rk,qj,qk,dqj,dqk
        real(real64), intent(out) :: value,dj(3),dk(3)
        real(real64) :: qjk,dqjk,rjk,ujk(3),cosine,a,da,product
        value=0;dj=0;dk=0
        if (rj <= 1e-12_real64 .or. rk <= 1e-12_real64 .or. rj > p(1) .or. rk > p(1)) return
        qjk=1;dqjk=0;ujk=0
        if (f(1) == 4) then
            ujk=rk*uk-rj*uj;rjk=sqrt(sum(ujk**2))
            if (rjk <= 1e-12_real64 .or. rjk > p(1)) return
            ujk=ujk/rjk
            call generic_radial(2,rjk,f(4),p,qjk,dqjk)
        end if
        cosine=max(-1.0_real64,min(1.0_real64,sum(uj*uk)))
        call angular_power(cosine,p(4),p(5),f(5),0.5_real64*p(5)*p(4),a,da)
        product=qj*qk*qjk
        value=2*a*product
        dj=2*(da*product*(uk-cosine*uj)/rj+a*qk*(dqj*qjk*uj-qj*dqjk*ujk))
        dk=2*(da*product*(uj-cosine*uk)/rk+a*qj*(dqk*qjk*uk+qk*dqjk*ujk))
    end subroutine
'''
s=s.replace('end module',helper+'end module');p.write_text(s)
p=root/'accelnet_target_kernels.f90';s=p.read_text()
s=s.replace('powers, radial_cache, use_moment','powers, radial_cache, jacobian, use_moment')
s=s.replace('powers(:,:,:), radial_cache(:,:,:)','powers(:,:,:), radial_cache(:,:,:), jacobian(:,:,:)')
s=s.replace('geom,radial_cache,g)', 'geom,radial_cache,jacobian,g)')
s=s.replace('geom,edge_row,radial_cache,g,edge_force)', 'geom,edge_row,radial_cache,jacobian,g,edge_force)')
start=s.index('    subroutine generic_values');stop=s.index('    subroutine generic_forces',start)
prefix=s[:start];values=s[start:stop];force=s[stop:]
values=values.replace('integer :: row,b,s,j,k,tj,tk,kind,t1,t2,group','integer :: row,b,s,j,k,tj,tk,kind,t1,t2,group,c')
values=values.replace('real(real64), contiguous, intent(inout) :: g(:,:)','real(real64), contiguous, intent(inout) :: g(:,:),jacobian(:,:,:)')
values=values.replace('gradient(3),total12','gradient(3),gradient_k(3),total12').replace('group,total,v,dv,gradient,total12','group,c,total,v,dv,gradient,gradient_k,total12')
values=values.replace('''                total = 0; total12 = 0
                do j''','''                total = 0; total12 = 0
                if (kind == 4 .or. kind == 5) then
                    do j = offsets(row),offsets(row+1)-1
                        jacobian(j,b,:) = 0
                    end do
                end if
                do j''')
values=values.replace('''                            v = generic_pair_value(features(:,b,s),fp(:,b,s),geom(1:3,j),geom(1:3,k), &
                                geom(4,j),geom(4,k),radial_cache(j,group,1),radial_cache(k,group,1))''','''                            call generic_pair_both(features(:,b,s),fp(:,b,s),geom(1:3,j),geom(1:3,k), &
                                geom(4,j),geom(4,k),radial_cache(j,group,1),radial_cache(k,group,1), &
                                radial_cache(j,group,2),radial_cache(k,group,2),v,gradient,gradient_k)
                            do c = 1,3
                                jacobian(j,b,c) = jacobian(j,b,c)+gradient(c)
                                jacobian(k,b,c) = jacobian(k,b,c)+gradient_k(c)
                            end do''')
force=force.replace('geom,edge_row,radial_cache,g,edge_force)', 'geom,edge_row,radial_cache,jacobian,g,edge_force)')
force=force.replace('fp(:,:,:),radial_cache(:,:,:),geom(:,:),g(:,:)', 'fp(:,:,:),radial_cache(:,:,:),jacobian(:,:,:),geom(:,:),g(:,:)')
start=force.index('                    if (tj /= t1 .and. tj /= t2) cycle');stop=force.index('                end if',start)
force=force[:start]+'''                    f = f+g(row,b)*jacobian(j,b,:)
'''+force[stop:]
p.write_text(prefix+values+force)
# Keep all source/directive lines within GNU's standard free-form limit.
for p in root.glob('*.f90'):
 out=[]
 for line in p.read_text().splitlines():
  while len(line)>130 and (not line.lstrip().startswith('!') or line.lstrip().startswith('!$omp')):
   k=line.rfind(',',0,110)
   if k<0:break
   out.append(line[:k+1]+' &')
   line=('        !$omp& ' if line.lstrip().startswith('!$omp') else '            ')+line[k+1:].lstrip()
  out.append(line)
 p.write_text('\n'.join(out)+'\n')
