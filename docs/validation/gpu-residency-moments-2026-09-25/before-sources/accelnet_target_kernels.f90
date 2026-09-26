! All arguments are plain contiguous arrays. No derived-type deep mapping,
! device allocation, or descriptor Jacobian is needed by the kernels.
module accelnet_target_kernels
    use iso_fortran_env, only: real64
    use omp_lib, only: omp_is_initial_device
    use accelnet_target_math
    implicit none
    private
    public :: run_target_batch
    real(real64), parameter :: pi = 3.14159265358979_real64, eps = 1e-12_real64
contains
    subroutine run_target_batch(device, meta, nodes, acts, woffset, weights, params, shift, scale, spin, &
                                species, centers, offsets, indices, dr, energies, forces, virial, g, values, deriv, delta)
        integer, intent(in) :: device
        integer, contiguous, intent(in) :: meta(:,:), nodes(:,:), acts(:,:), woffset(:,:)
        real(real64), contiguous, intent(in) :: weights(:,:), params(:,:), shift(:,:), scale(:,:), spin(:,:)
        integer, contiguous, intent(in) :: species(:), centers(:), offsets(:), indices(:)
        real(real64), contiguous, intent(in) :: dr(:,:)
        real(real64), contiguous, intent(out) :: energies(:)
        real(real64), contiguous, intent(inout) :: forces(:,:), virial(:,:)
        real(real64), contiguous, intent(inout) :: g(:,:), values(:,:,:), deriv(:,:,:), delta(:,:,:)
        integer :: row, s, nr, na, dim, multi, version, ct, j, k, b, l, i, o, nin, nout, a, c, target, nrw
        real(real64) :: rj, rk, fcj, fck, dfcj, sj, sk, cosine, x, t0, t1, t2, dt0, dt1, dt2
        real(real64) :: v, dv, z, rc, ac, alpha, xscale, vc, dc, coeff, radial, fj(3), center_f(3), w(3,3)
        logical :: on_host
        nrw = size(centers)
        on_host = .true.
        ! One structured data lifetime spans descriptor, NN/backprop, contraction
        ! and scatter. Only inputs and final results cross the host boundary.
        !$omp target data device(device) &
        !$omp& map(to: meta, nodes, acts, woffset, weights, params, shift, scale, spin) &
        !$omp& map(to: species, centers, offsets, indices, dr) &
        !$omp& map(from: energies) map(tofrom: forces, virial) map(alloc: g, values, deriv, delta)
        !$omp target device(device) map(from: on_host)
        on_host = omp_is_initial_device()
        !$omp end target
        if (on_host) error stop 'OpenMP target: CPU fallback is not GPU execution'

        !$omp target teams distribute parallel do device(device) &
        !$omp& private(s,nr,na,dim,multi,version,ct,j,k,b,l,i,o,nin,nout,rj,rk,fcj,fck,sj,sk,cosine,x, &
        !$omp& t0,t1,t2,v,z,rc,ac,alpha,xscale)
        do row = 1, nrw
            s = species(centers(row))
            nr = meta(1,s); na = meta(2,s); multi = meta(3,s); version = meta(4,s); ct = meta(5,s)
            dim = nodes(1,s); rc = params(1,s); ac = params(2,s); alpha = params(3,s)
            do b = 1, dim
                g(row,b) = 0.0_real64
            end do
            if (version == 10) then
                sj = 0.0_real64
                if (multi == 1) sj = spin(species(indices(offsets(row)+meta(7,s)-1)),s)
                do b = 1, nr
                    v = real((-1)**(b-1),real64)
                    g(row,b) = v
                    if (multi == 1) g(row,nr+na+b) = sj*v
                end do
            end if
            do j = offsets(row), offsets(row+1)-1
                rj = sqrt(sum(dr(:,j)**2))
                sj = spin(species(indices(j)),s)
                if (rj <= rc .and. rj > eps) then
                    fcj = target_cutoff_value(rj,rc,ct,alpha)
                    x = 2.0_real64*rj/rc-1.0_real64
                    t0 = 1.0_real64; t1 = x
                    do b = 1, nr
                        v = fcj*t0
                        g(row,b) = g(row,b)+v
                        if (multi == 1) g(row,nr+na+b) = g(row,nr+na+b)+sj*v
                        t2 = 2.0_real64*x*t1-t0; t0 = t1; t1 = t2
                    end do
                end if
                if (rj > ac .or. rj < eps) cycle
                fcj = target_cutoff_value(rj,ac,ct,alpha)
                do k = j+1, offsets(row+1)-1
                    rk = sqrt(sum(dr(:,k)**2))
                    if (rk > ac .or. rk < eps) cycle
                    fck = target_cutoff_value(rk,ac,ct,alpha)
                    sk = spin(species(indices(k)),s)
                    cosine = sum(dr(:,j)*dr(:,k))/(rj*rk)
                    x = cosine
                    if (version == 1) x = (2.0_real64*cosine-pi)/pi
                    t0 = 1.0_real64; t1 = x
                    do b = 1, na
                        v = fcj*fck*t0
                        g(row,nr+b) = g(row,nr+b)+v
                        if (multi == 1) g(row,2*nr+na+b) = g(row,2*nr+na+b)+sj*sk*v
                        t2 = 2.0_real64*x*t1-t0; t0 = t1; t1 = t2
                    end do
                end do
            end do
            ! NN work is laid out with batch rows contiguous, including networks
            ! with different depths and widths for different elements.
            do i = 1, dim
                values(row,i,1) = (g(row,i)-shift(i,s))*scale(i,s)
            end do
            do l = 1, meta(6,s)-1
                nin = nodes(l,s); nout = nodes(l+1,s); o = woffset(l,s)+1
                do j = 1, nout
                    z = weights(o+nin*nout+j-1,s)
                    do i = 1, nin
                        z = z+weights(o+(i-1)*nout+j-1,s)*values(row,i,l)
                    end do
                    v = target_activate(z,acts(l,s))
                    values(row,j,l+1) = v
                    deriv(row,j,l+1) = target_activation_derivative(z,v,acts(l,s))
                end do
            end do
            energies(row) = values(row,1,meta(6,s))/params(4,s)+params(5,s)
            delta(row,1,1) = 1.0_real64
            do l = meta(6,s)-1, 1, -1
                nin = nodes(l,s); nout = nodes(l+1,s); o = woffset(l,s)+1
                do j = 1, nout
                    delta(row,j,1) = delta(row,j,1)*deriv(row,j,l+1)
                end do
                do i = 1, nin
                    v = 0.0_real64
                    do j = 1, nout
                        v = v+weights(o+(i-1)*nout+j-1,s)*delta(row,j,1)
                    end do
                    delta(row,i,2) = v
                end do
                do i = 1, nin
                    delta(row,i,1) = delta(row,i,2)
                end do
            end do
            do b = 1, dim
                g(row,b) = -delta(row,b,1)*scale(b,s)/params(4,s)
            end do
        end do
        !$omp end target teams distribute parallel do

        !$omp target teams distribute parallel do device(device) &
        !$omp& private(s,nr,na,multi,version,ct,j,k,b,a,c,target,rj,rk,fcj,fck,dfcj,sj,sk,cosine,x, &
        !$omp& t0,t1,t2,dt0,dt1,dt2,v,dv,rc,ac,alpha,xscale,vc,dc,coeff,radial,fj,center_f,w)
        do row = 1, nrw
            s = species(centers(row))
            nr = meta(1,s); na = meta(2,s); multi = meta(3,s); version = meta(4,s); ct = meta(5,s)
            rc = params(1,s); ac = params(2,s); alpha = params(3,s)
            center_f = 0.0_real64; w = 0.0_real64
            ! Revisit geometry and contract immediately. For each edge j, sum
            ! all angular partners k before one atomic scatter per component.
            ! There is no (atom,descriptor,neighbor,xyz) Jacobian allocation.
            do j = offsets(row), offsets(row+1)-1
                rj = sqrt(sum(dr(:,j)**2)); sj = spin(species(indices(j)),s)
                fj = 0.0_real64
                if (rj <= rc .and. rj > eps) then
                    fcj = target_cutoff_value(rj,rc,ct,alpha)
                    dfcj = target_cutoff_derivative(rj,rc,ct,alpha)
                    x = 2.0_real64*rj/rc-1.0_real64; xscale = 2.0_real64/rc
                    t0 = 1.0_real64; t1 = x; dt0 = 0.0_real64; dt1 = xscale; radial = 0.0_real64
                    do b = 1, nr
                        coeff = g(row,b)
                        if (multi == 1) coeff = coeff+sj*g(row,nr+na+b)
                        radial = radial+coeff*(dfcj*t0+fcj*dt0)
                        t2 = 2*x*t1-t0; dt2 = 2*xscale*t1+2*x*dt1-dt0
                        t0 = t1; t1 = t2; dt0 = dt1; dt1 = dt2
                    end do
                    fj = radial*dr(:,j)/rj
                end if
                if (rj <= ac .and. rj >= eps) then
                    fcj = target_cutoff_value(rj,ac,ct,alpha)
                    dfcj = target_cutoff_derivative(rj,ac,ct,alpha)
                    do k = offsets(row), offsets(row+1)-1
                        if (k == j) cycle
                        rk = sqrt(sum(dr(:,k)**2))
                        if (rk > ac .or. rk < eps) cycle
                        fck = target_cutoff_value(rk,ac,ct,alpha); sk = spin(species(indices(k)),s)
                        cosine = sum(dr(:,j)*dr(:,k))/(rj*rk)
                        x = cosine; xscale = 1.0_real64
                        if (version == 1) then
                            x = (2*cosine-pi)/pi; xscale = 2.0_real64/pi
                        end if
                        t0 = 1.0_real64; t1 = x; dt0 = 0.0_real64; dt1 = xscale
                        vc = 0.0_real64; dc = 0.0_real64
                        do b = 1, na
                            coeff = g(row,nr+b)
                            if (multi == 1) coeff = coeff+sj*sk*g(row,2*nr+na+b)
                            vc = vc+coeff*t0; dc = dc+coeff*dt0
                            t2 = 2*x*t1-t0; dt2 = 2*xscale*t1+2*x*dt1-dt0
                            t0 = t1; t1 = t2; dt0 = dt1; dt1 = dt2
                        end do
                        fj = fj+dfcj*fck*vc*dr(:,j)/rj + &
                            fcj*fck*dc*(dr(:,k)/(rj*rk)-cosine*dr(:,j)/(rj*rj))
                    end do
                end if
                target = indices(j)
                do c = 1, 3
                    !$omp atomic update
                    forces(c,target) = forces(c,target)+fj(c)
                    center_f(c) = center_f(c)-fj(c)
                    do a = 1, 3
                        w(a,c) = w(a,c)+dr(a,j)*fj(c)
                    end do
                end do
            end do
            target = centers(row)
            do c = 1, 3
                !$omp atomic update
                forces(c,target) = forces(c,target)+center_f(c)
                do a = 1, 3
                    !$omp atomic update
                    virial(a,c) = virial(a,c)+w(a,c)
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
        !$omp end target data
    end subroutine
end module
