! All arguments are plain contiguous arrays. No derived-type deep mapping,
! device allocation, or descriptor Jacobian is needed by the kernels.
module accelnet_target_kernels
    use iso_fortran_env, only: real64
    use omp_lib, only: omp_get_wtime
    use accelnet_target_math
    implicit none
    private
    public :: run_target_batch
    real(real64), parameter :: pi = 3.14159265358979_real64, eps = 1e-12_real64
contains
    subroutine run_target_batch(device, meta, nodes, acts, woffset, weights, params, shift, scale, spin, &
                                mp, multiplicity, polynomial, species, centers, offsets, indices, dr, &
                                energies, forces, virial, g, values, deriv, delta, moments, powers, use_moment, &
                                geom, edge_row, edge_force, nrw, natoms, nedges, stages)
        integer, intent(in) :: device, nrw, natoms, nedges
        integer, contiguous, intent(in) :: meta(:,:), nodes(:,:), acts(:,:), woffset(:,:), mp(:,:,:)
        real(real64), contiguous, intent(in) :: multiplicity(:,:), polynomial(:,:,:)
        real(real64), contiguous, intent(in) :: weights(:,:), params(:,:), shift(:,:), scale(:,:), spin(:,:)
        integer, contiguous, intent(in) :: species(:), centers(:), offsets(:), indices(:)
        real(real64), contiguous, intent(in) :: dr(:,:)
        real(real64), contiguous, intent(out) :: energies(:)
        real(real64), contiguous, intent(inout) :: forces(:,:), virial(:,:)
        real(real64), contiguous, intent(inout) :: g(:,:), values(:,:,:), deriv(:,:,:), delta(:,:,:)
        real(real64), contiguous, intent(inout) :: moments(:,:,:), powers(:,:,:)
        integer, contiguous, intent(inout) :: use_moment(:), edge_row(:)
        real(real64), contiguous, intent(inout) :: geom(:,:), edge_force(:,:)
        real(real64), intent(out) :: stages(3)
        integer :: row, s, nr, na, dim, multi, version, ct, j, k, b, l, i, o, nin, nout, a, c, target
        integer :: entry, q, ax, ay, az, angular_neighbors, nm, tile, lane
        real(real64) :: rj, rk, fcj, fck, dfcj, sj, sk, cosine, x, t0, t1, t2, dt0, dt1, dt2
        real(real64) :: v, dv, z, rc, ac, alpha, xscale, vc, dc, coeff, radial, fj(3), center_f(3), w(3,3)
        real(real64) :: self0, self1, mono, u(3), grad(3), dm(3), correction, started
        ! All arrays have persistent mappings owned by target_workspace.
        ! Structured references here neither copy data nor allocate device memory.
        !$omp target data device(device) &
        !$omp& map(alloc: meta, nodes, acts, woffset, weights, params, shift, scale, spin, mp, multiplicity, polynomial) &
        !$omp& map(alloc: species, centers, offsets, indices, dr, energies, forces, virial, &
        !$omp& g, values, deriv, delta, moments, powers, use_moment, geom, edge_row, edge_force)
        started = omp_get_wtime()
        !$omp target teams distribute parallel do device(device) thread_limit(32) &
        !$omp& private(s,nr,na,dim,multi,version,ct,j,k,b,l,i,o,nin,nout,rj,rk,fcj,fck,sj,sk,cosine,x, &
        !$omp& t0,t1,t2,v,z,rc,ac,alpha,xscale,entry,q,ax,ay,az,nm,angular_neighbors,self0,self1,mono,u)
        do row = 1, max(nrw,natoms)
            if (row <= natoms) forces(:,row) = 0.0_real64
            if (row == 1) virial = 0.0_real64
            if (row > nrw) cycle
            s = species(centers(row))
            nr = meta(1,s); na = meta(2,s); multi = meta(3,s); version = meta(4,s); ct = meta(5,s)
            dim = nodes(1,s); rc = params(1,s); ac = params(2,s); alpha = params(3,s)
            angular_neighbors = 0
            do j = offsets(row), offsets(row+1)-1
                rj = sqrt(sum(dr(:,j)**2))
                if (rj <= ac .and. rj > eps) angular_neighbors = angular_neighbors+1
                edge_row(j) = row
                geom(1:3,j) = 0
                if (rj >= eps) geom(1:3,j) = dr(:,j)/rj
                geom(4,j) = rj
                geom(5,j) = target_cutoff_value(rj,ac,ct,alpha)
                geom(6,j) = target_cutoff_derivative(rj,ac,ct,alpha)
                geom(7,j) = spin(species(indices(j)),s)
            end do
            use_moment(row) = 0
            ! Initial operation-count heuristic. Forced direct/moment modes are
            ! available for measurement; do not inherit the CPU's threshold.
            if (meta(8,s) == 2 .or. (meta(8,s) == 0 .and. &
                real(angular_neighbors,real64)*na >= real(meta(9,s),real64))) use_moment(row) = 1
            nm = meta(9,s)
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
                rj = geom(4,j)
                sj = geom(7,j)
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
                fcj = geom(5,j)
                if (use_moment(row) == 1) then
                    u = geom(1:3,j)
                    powers(j,1,:) = 1
                    do q = 2, na
                        powers(j,q,:) = powers(j,q-1,:)*u
                    end do
                    cycle
                end if
                do k = j+1, offsets(row+1)-1
                    rk = geom(4,k)
                    if (rk > ac .or. rk < eps) cycle
                    fck = geom(5,k); sk = geom(7,k)
                    cosine = sum(geom(1:3,j)*geom(1:3,k))
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
        end do
        !$omp end target teams distribute parallel do
        if (any(meta(8,:) /= 1)) then
            ! Parallelize over both rows and monomials, instead of leaving one
            ! GPU thread to build every moment for an atom.
            !$omp target teams distribute parallel do collapse(3) device(device) &
            !$omp& private(row,s,j,ax,ay,az,mono,fcj,sj,self0,self1)
            do tile = 0, (nrw+512-1)/512-1
              do entry = 1, size(mp,2)
                do lane = 1, 512
                    row = tile*512+lane
                    if (row > nrw) cycle
                    if (use_moment(row) /= 1) cycle
                    s = species(centers(row))
                    if (entry > meta(9,s)) cycle
                    ax = mp(1,entry,s); ay = mp(2,entry,s); az = mp(3,entry,s)
                    self0 = 0; self1 = 0
                    do j = offsets(row), offsets(row+1)-1
                        if (geom(4,j) > params(2,s) .or. geom(4,j) < eps) cycle
                        fcj = geom(5,j); sj = geom(7,j)
                        mono = powers(j,ax+1,1)*powers(j,ay+1,2)*powers(j,az+1,3)
                        self0 = self0+fcj*mono; self1 = self1+sj*fcj*mono
                    end do
                    moments(row,entry,1) = self0; moments(row,entry,2) = self1
                end do
            end do
            end do
            !$omp end target teams distribute parallel do
            !$omp target teams distribute parallel do device(device) thread_limit(32) &
            !$omp& private(s,nr,na,multi,nm,entry,q,b,j,self0,self1,fcj,sj)
            do row = 1, nrw
                if (use_moment(row) /= 1) cycle
                s = species(centers(row)); nr = meta(1,s); na = meta(2,s); multi = meta(3,s); nm = meta(9,s)
                self0 = 0; self1 = 0
                do j = offsets(row), offsets(row+1)-1
                    if (geom(4,j) > params(2,s) .or. geom(4,j) < eps) cycle
                    fcj = geom(5,j); sj = geom(7,j)
                    self0 = self0+fcj*fcj; self1 = self1+(sj*fcj)**2
                end do
                do q = 1, na
                    delta(row,q,1) = 0; delta(row,q,2) = 0
                end do
                do entry = 1, nm
                    q = mp(4,entry,s)
                    delta(row,q,1) = delta(row,q,1)+multiplicity(entry,s)*moments(row,entry,1)**2
                    delta(row,q,2) = delta(row,q,2)+multiplicity(entry,s)*moments(row,entry,2)**2
                end do
                do q = 1, na
                    delta(row,q,1) = 0.5_real64*(delta(row,q,1)-self0)
                    delta(row,q,2) = 0.5_real64*(delta(row,q,2)-self1)
                end do
                do b = 1, na
                    do q = 1, b
                        g(row,nr+b) = g(row,nr+b)+polynomial(b,q,s)*delta(row,q,1)
                        if (multi == 1) g(row,2*nr+na+b) = g(row,2*nr+na+b)+polynomial(b,q,s)*delta(row,q,2)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
        end if
        stages(1) = omp_get_wtime()-started
        started = omp_get_wtime()
        !$omp target teams distribute parallel do device(device) thread_limit(32) &
        !$omp& private(s,dim,l,nin,nout,o,j,i,z,v,b,q,nr,na,multi)
        do row = 1, nrw
            s = species(centers(row)); dim = nodes(1,s)
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
            nr = meta(1,s); na = meta(2,s); multi = meta(3,s)
            if (use_moment(row) == 1) then
                do q = 1, na
                    delta(row,q,1) = 0; delta(row,q,2) = 0
                    do b = q, na
                        delta(row,q,1) = delta(row,q,1)+g(row,nr+b)*polynomial(b,q,s)
                        if (multi == 1) delta(row,q,2) = delta(row,q,2)+g(row,2*nr+na+b)*polynomial(b,q,s)
                    end do
                end do
            end if
        end do
        !$omp end target teams distribute parallel do
        stages(2) = omp_get_wtime()-started
        started = omp_get_wtime()

        !$omp target teams distribute parallel do device(device) &
        !$omp& private(row,s,nr,na,multi,version,ct,k,b,rj,rk,fcj,fck,dfcj,sj,sk,cosine,x, &
        !$omp& t0,t1,t2,dt0,dt1,dt2,rc,ac,alpha,xscale,vc,dc,coeff,radial,fj, &
        !$omp& entry,q,ax,ay,az,nm,mono,u,grad,dm,correction)
        do j = 1, nedges
            row = edge_row(j); s = species(centers(row))
            nr = meta(1,s); na = meta(2,s); multi = meta(3,s); version = meta(4,s); ct = meta(5,s)
            rc = params(1,s); ac = params(2,s); alpha = params(3,s); nm = meta(9,s)
            ! Each edge contracts independently; a later row kernel scatters.
            ! Cache geometry and powers; never materialize a descriptor Jacobian.
                rj = geom(4,j); sj = geom(7,j)
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
                    fcj = geom(5,j); dfcj = geom(6,j)
                    if (use_moment(row) == 1) then
                        u = geom(1:3,j)
                        do entry = 1, nm
                            ax = mp(1,entry,s); ay = mp(2,entry,s); az = mp(3,entry,s); q = mp(4,entry,s)
                            mono = powers(j,ax+1,1)*powers(j,ay+1,2)*powers(j,az+1,3)
                            grad = 0
                            if (ax > 0) grad(1) = ax*powers(j,ax,1)*powers(j,ay+1,2)*powers(j,az+1,3)
                            if (ay > 0) grad(2) = ay*powers(j,ax+1,1)*powers(j,ay,2)*powers(j,az+1,3)
                            if (az > 0) grad(3) = az*powers(j,ax+1,1)*powers(j,ay+1,2)*powers(j,az,3)
                            dm = dfcj*u*mono+fcj*(grad-u*sum(u*grad))/rj
                            coeff = multiplicity(entry,s)*(delta(row,q,1)*moments(row,entry,1) + &
                                sj*delta(row,q,2)*moments(row,entry,2))
                            fj = fj+coeff*dm
                            ! Entries are ordered by degree; apply the self-image
                            ! correction after the last monomial of each degree.
                            if (ax == q-1 .and. ay == 0 .and. az == 0) then
                                correction = fcj*dfcj*(delta(row,q,1)+sj*sj*delta(row,q,2))
                                fj = fj-correction*u
                            end if
                        end do
                    else
                    do k = offsets(row), offsets(row+1)-1
                        if (k == j) cycle
                        rk = geom(4,k)
                        if (rk > ac .or. rk < eps) cycle
                        fck = geom(5,k); sk = geom(7,k)
                        cosine = sum(geom(1:3,j)*geom(1:3,k))
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
                end if
            edge_force(:,j) = fj
        end do
        !$omp end target teams distribute parallel do
        !$omp target teams distribute parallel do device(device) thread_limit(32) &
        !$omp& private(j,target,c,a,fj,center_f,w)
        do row = 1, nrw
            center_f = 0; w = 0
            do j = offsets(row), offsets(row+1)-1
                fj = edge_force(:,j)
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
        stages(3) = omp_get_wtime()-started
        !$omp end target data
    end subroutine
end module
