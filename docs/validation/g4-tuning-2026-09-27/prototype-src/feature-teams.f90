! All arguments are plain contiguous arrays. No derived-type deep mapping,
! device allocation is performed inside the kernels. G4 saves its Jacobian.
module accelnet_target_kernels
    use iso_fortran_env, only: real64
    use accelnet_target_runtime, only: omp_get_wtime, omp_get_initial_device
    use accelnet_target_math
    implicit none
    private
    public :: run_target_batch
    real(real64), parameter :: pi = 3.14159265358979_real64, eps = 1e-12_real64
contains
    subroutine run_target_batch(device, meta, nodes, acts, woffset, weights, params, shift, scale, spin, &
                                features, feature_params, local_species, mp, multiplicity, polynomial, &
                                species, centers, offsets, indices, dr, &
                                energies, forces, virial, g, values, deriv, delta, moments, powers, &
            radial_cache, jacobian, g4_scratch, g4_first, use_moment, &
                                geom, edge_row, edge_force, nrw, natoms, nedges, stages)
        integer, contiguous, intent(in) :: features(:,:,:),local_species(:,:)
        real(real64), contiguous, intent(in) :: feature_params(:,:,:)
        integer, intent(in) :: device, nrw, natoms, nedges
        integer, contiguous, intent(in) :: meta(:,:), nodes(:,:), acts(:,:), woffset(:,:), mp(:,:,:)
        real(real64), contiguous, intent(in) :: multiplicity(:,:), polynomial(:,:,:)
        real(real64), contiguous, intent(in) :: weights(:,:), params(:,:), shift(:,:), scale(:,:), spin(:,:)
        integer, contiguous, intent(in) :: species(:), centers(:), offsets(:), indices(:)
        real(real64), contiguous, intent(in) :: dr(:,:)
        real(real64), contiguous, intent(out) :: energies(:)
        real(real64), contiguous, intent(inout) :: forces(:,:), virial(:,:)
        real(real64), contiguous, intent(inout) :: g(:,:), values(:,:,:), deriv(:,:,:), delta(:,:,:)
        real(real64), contiguous, intent(inout) :: moments(:,:,:), powers(:,:,:), radial_cache(:,:,:), &
            jacobian(:,:,:), g4_scratch(:,:,:)
        integer, contiguous, intent(inout) :: use_moment(:), edge_row(:), g4_first(:,:,:)
        real(real64), contiguous, intent(inout) :: geom(:,:), edge_force(:,:)
        real(real64), intent(out) :: stages(3)
        integer :: row, s, nr, na, dim, multi, version, ct, j, k, b, l, i, o, nin, nout, a, c, target
        integer :: entry, q, ax, ay, az, angular_neighbors, nm, bb
        real(real64) :: rj, rk, fcj, fck, dfcj, sj, sk, cosine, x, t0, t1, t2, dt0, dt1, dt2
        real(real64) :: v, dv, z, rc, ac, alpha, xscale, vc, dc, coeff, radial, fj(3), center_f(3), w(3,3)
        real(real64) :: self0, self1, mono, u(3), grad(3), dm(3), correction, started
        real(real64) :: hy, hz, dhy, dhz, hyz
        ! All arrays have persistent mappings owned by target_workspace.
        ! Structured references here neither copy data nor allocate device memory.
        !$omp target data device(device) if(device /= omp_get_initial_device()) &
        !$omp& map(alloc: features,feature_params,local_species) &
        !$omp& map(alloc: meta, nodes, acts, woffset, weights, params, shift, scale, spin, mp, multiplicity, polynomial) &
        !$omp& map(alloc: species, centers, offsets, indices, dr, energies, forces, virial, &
        !$omp& g, values, deriv, delta, moments, powers, radial_cache, jacobian, g4_scratch, g4_first, &
        !$omp& use_moment, geom, edge_row, edge_force)
        started = omp_get_wtime()
        !$omp target teams distribute parallel do device(device) if(device /= omp_get_initial_device()) thread_limit(32) &
        !$omp& private(s,nr,na,dim,multi,version,ct,j,k,b,l,i,o,nin,nout,rj,rk,fcj,fck,sj,sk,cosine,x, &
        !$omp& t0,t1,t2,v,z,rc,ac,alpha,xscale,entry,q,ax,ay,az,nm,angular_neighbors,self0,self1,mono,u)
        do row = 1, max(nrw,natoms)
            if (row <= natoms) forces(:,row) = 0.0_real64
            if (row == 1) virial = 0.0_real64
            if (row > nrw) cycle
            s = species(centers(row))
            nr = meta(1,s); na = meta(2,s); multi = meta(3,s); version = meta(4,s); ct = meta(5,s)
            dim = nodes(1,s); rc = params(1,s); ac = params(2,s); alpha = params(3,s)
            if (meta(10,s) == 1) then
                use_moment(row) = 0
                do j = offsets(row), offsets(row+1)-1
                    edge_row(j) = row
                    rj = sqrt(sum(dr(:,j)**2)); geom(4,j) = rj
                    geom(1:3,j) = 0
                    if (rj > eps) geom(1:3,j) = dr(:,j)/rj
                end do
                ! Cache q and q' once per radial parameter group and edge.
                ! Fused with geometry to avoid another GPU launch. Chebyshev
                ! and generic rows own disjoint edges in mixed-family models.
                do b = 1,nodes(1,s)
                    if (features(8,b,s) /= b) cycle
                    o = features(7,b,s)
                    do j = offsets(row),offsets(row+1)-1
                        call generic_radial(2,geom(4,j),features(4,b,s),feature_params(:,b,s), &
                            radial_cache(j,o,1),radial_cache(j,o,2))
                    end do
                end do
                cycle
            end if
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
            ! Keep one atom's monomials adjacent in the GPU work order so their
            ! neighbor geometry and powers can be reused before cache eviction.
            !$omp target teams distribute parallel do collapse(2) device(device) if(device /= omp_get_initial_device()) &
            !$omp& private(s,j,ax,ay,az,mono,fcj,sj,self0,self1)
            do row = 1, nrw
                do entry = 1, size(mp,2)
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
            !$omp end target teams distribute parallel do
            !$omp target teams distribute parallel do device(device) if(device /= omp_get_initial_device()) thread_limit(32) &
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
        if (any(features(1,:,:) == 4)) call g4_values_derivatives(device,nrw,meta,nodes,features,feature_params, &
            local_species,species,centers,offsets,indices,dr,geom,radial_cache,g,jacobian,g4_scratch,g4_first)
        if (any(features(1,:,:) > 0 .and. features(1,:,:) /= 4 .and. features(1,:,:) /= 7)) then
            if (any(features(12,:,:) > 0 .and. features(1,:,:) == 5)) then
                call generic_values_grouped(device,nrw,meta,nodes,features,feature_params, &
                    local_species,species,centers,offsets,indices,geom,radial_cache,g)
            else
                call generic_values(device,nrw,meta,nodes,features,feature_params, &
                    local_species,species,centers,offsets,indices,geom,radial_cache,g)
            end if
        end if
        stages(1) = omp_get_wtime()-started
        started = omp_get_wtime()
        !$omp target teams distribute parallel do device(device) if(device /= omp_get_initial_device()) thread_limit(32) &
        !$omp& private(s,dim,l,nin,nout,o,j,i,z,v,b,q,nr,na,multi,entry,self0,self1,bb)
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
            ! Reuse NN gradients for equal-power coefficient sums. Distinct
            ! powers already have the desired coefficient, so leave them in place.
            ! No additional array, allocation, transfer, or kernel launch.
            if (meta(10,s) == 1) then
                do b=1,dim
                    if (features(1,b,s) == 4) cycle ! G4 contracts its saved Jacobian after the NN.
                    if (features(13,b,s) /= b .or. features(14,b,s) == 0) cycle
                    v=g(row,b); bb=features(14,b,s)
                    do while (bb /= 0)
                        v=v+g(row,bb)
                        bb=features(14,bb,s)
                    end do
                    g(row,b)=v
                end do
            end if
            nr = meta(1,s); na = meta(2,s); multi = meta(3,s)
            if (use_moment(row) == 1) then
                do q = 1, na
                    delta(row,q,1) = 0; delta(row,q,2) = 0
                    do b = q, na
                        delta(row,q,1) = delta(row,q,1)+g(row,nr+b)*polynomial(b,q,s)
                        if (multi == 1) delta(row,q,2) = delta(row,q,2)+g(row,2*nr+na+b)*polynomial(b,q,s)
                    end do
                end do
                ! Raw moments are dead after descriptor evaluation. Reuse their
                ! storage for the two species-weight channels of the force
                ! polynomial, once per center rather than once per edge.
                do entry = 1, meta(9,s)
                    q = mp(4,entry,s)
                    moments(row,entry,1) = multiplicity(entry,s)*delta(row,q,1)*moments(row,entry,1)
                    moments(row,entry,2) = multiplicity(entry,s)*delta(row,q,2)*moments(row,entry,2)
                end do
                self0 = 0; self1 = 0
                do q = 1, na
                    self0 = self0+delta(row,q,1)
                    self1 = self1+delta(row,q,2)
                end do
                ! The self-pair correction is independent of direction.
                delta(row,1,1) = self0; delta(row,1,2) = self1
            end if
        end do
        !$omp end target teams distribute parallel do
        stages(2) = omp_get_wtime()-started
        started = omp_get_wtime()

        !$omp target teams distribute parallel do device(device) if(device /= omp_get_initial_device()) &
        !$omp& private(row,s,nr,na,multi,version,ct,k,b,rj,rk,fcj,fck,dfcj,sj,sk,cosine,x, &
        !$omp& t0,t1,t2,dt0,dt1,dt2,rc,ac,alpha,xscale,vc,dc,coeff,radial,fj, &
        !$omp& entry,q,ax,ay,az,nm,mono,u,grad,dm,correction,v,hy,hz,dhy,dhz,hyz)
        do j = 1, nedges
            row = edge_row(j); s = species(centers(row))
            if (meta(10,s) == 1) cycle
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
                    fj = radial*geom(1:3,j)
                end if
                if (rj <= ac .and. rj >= eps) then
                    fcj = geom(5,j); dfcj = geom(6,j)
                    if (use_moment(row) == 1) then
                        u = geom(1:3,j)
                        ! Evaluate the total-degree polynomial and all three
                        ! partial derivatives together. Differentiating Horner
                        ! uses the previous value before each multiply/add and
                        ! remains valid when any component of u is zero.
                        v = 0; dm = 0; q = nm
                        do ax = na-1, 0, -1
                            hy = 0; dhy = 0; hyz = 0
                            do ay = na-1-ax, 0, -1
                                hz = 0; dhz = 0
                                do az = na-1-ax-ay, 0, -1
                                    entry = q; q = q-1
                                    coeff = moments(row,entry,1)+sj*moments(row,entry,2)
                                    dhz = dhz*u(3)+hz
                                    hz = hz*u(3)+coeff
                                end do
                                dhy = dhy*u(2)+hy
                                hyz = hyz*u(2)+dhz
                                hy = hy*u(2)+hz
                            end do
                            dm(1) = dm(1)*u(1)+v
                            dm(2) = dm(2)*u(1)+dhy
                            dm(3) = dm(3)*u(1)+hyz
                            v = v*u(1)+hy
                        end do
                        ! Contract the polynomial and its Cartesian gradient
                        ! before projecting onto the tangent plane of u.
                        correction = fcj*dfcj*(delta(row,1,1)+sj*sj*delta(row,1,2))
                        fj = fj+(dfcj*v-correction)*u+fcj*(dm-u*sum(u*dm))/rj
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
                        ! Clenshaw evaluates the contracted series and its
                        ! derivative without forming every T_b and T'_b.
                        t0 = 0; t1 = 0; dt0 = 0; dt1 = 0
                        do b = na, 2, -1
                            coeff = g(row,nr+b)
                            if (multi == 1) coeff = coeff+sj*sk*g(row,2*nr+na+b)
                            t2 = coeff+2*x*t0-t1
                            dt2 = 2*t0+2*x*dt0-dt1
                            t1 = t0; t0 = t2; dt1 = dt0; dt0 = dt2
                        end do
                        coeff = g(row,nr+1)
                        if (multi == 1) coeff = coeff+sj*sk*g(row,2*nr+na+1)
                        vc = coeff+x*t0-t1
                        dc = xscale*(t0+x*dt0-dt1)
                        ! Geometry already contains unit directions. Reuse them
                        ! instead of normalizing both vectors for every pair.
                        fj = fj+(dfcj*fck*vc)*geom(1:3,j) + &
                            (fcj*fck*dc/rj)*(geom(1:3,k)-cosine*geom(1:3,j))
                    end do
                    end if
                end if
            edge_force(:,j) = fj
        end do
        !$omp end target teams distribute parallel do
        if (any(meta(10,:) == 1)) call generic_forces(device,nedges,meta,nodes,features,feature_params, &
            local_species,species,centers,offsets,indices,geom,edge_row,radial_cache,jacobian,g,edge_force)
        !$omp target teams distribute parallel do device(device) if(device /= omp_get_initial_device()) thread_limit(32) &
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
    ! The same center -> unordered pair -> descriptor traversal as the CPU
    ! G4 evaluator: values and both derivatives are accumulated together. Each
    ! center owns its CSR edges, so neither CPU nor GPU needs Jacobian atomics.
    subroutine g4_values_derivatives(device,nrw,meta,nodes,features,fp,local_species,species,centers,offsets, &
            indices,dr,geom,radial_cache,g,jacobian,g4_scratch,g4_first)
        integer, intent(in) :: device,nrw
        integer, contiguous, intent(in) :: meta(:,:),nodes(:,:),features(:,:,:),local_species(:,:)
        integer, contiguous, intent(in) :: species(:),centers(:),offsets(:),indices(:)
        real(real64), contiguous, intent(in) :: fp(:,:,:),dr(:,:),geom(:,:),radial_cache(:,:,:)
        real(real64), contiguous, intent(inout) :: g(:,:),jacobian(:,:,:),g4_scratch(:,:,:)
        integer, contiguous, intent(inout) :: g4_first(:,:,:)
        integer :: row,s,dim,b,j,k,tj,tk,t1,t2,head,group,cr,er,ar,c,last_group,last_cr,last_er,last_ar,q
        real(real64) :: rcmax,rj,rk,rjk,rjk2,cosine,uj(3),uk(3),ujk(3),dcj(3),dck(3)
        real(real64) :: qj,qk,qjk,dqj,dqk,dqjk,a,da,ca,cv,dj(3),dk(3),fc,dfc,ex,radial,rgradj(3),rgradk(3)
        !$omp target teams distribute device(device) if(device /= omp_get_initial_device()) thread_limit(32) num_teams(max(1,nrw)) &
        !$omp& map(alloc:meta,nodes,features,fp,local_species,species,centers,offsets,indices, &
        !$omp& dr,geom,radial_cache,g,jacobian,g4_scratch,g4_first) &
        !$omp& private(s,dim,b,j,k,tj,tk,t1,t2,head,group,cr,er,ar,c,last_group,last_cr,last_er,last_ar,q,rcmax,rj,rk,rjk,rjk2, &
        !$omp& cosine,uj,uk,ujk,dcj,dck,qj,qk,qjk,dqj,dqk,dqjk,a,da,ca,cv,dj,dk,fc,dfc,ex,radial,rgradj,rgradk)
        do row = 1,nrw
            s = species(centers(row))
            if (meta(10,s) /= 1) cycle
            dim = nodes(1,s); g4_first(:,:,row) = 0; rcmax = 0
            do b = dim,1,-1
                if (features(1,b,s) /= 4) cycle
                t1 = features(2,b,s); t2 = features(3,b,s)
                g4_first(t1,t2,row) = b; g4_first(t2,t1,row) = b
                rcmax = max(rcmax,fp(1,b,s))
                g(row,b) = 0
                do j = offsets(row),offsets(row+1)-1
                    jacobian(:,b,j) = 0
                end do
            end do
            if (rcmax == 0) cycle
            !$omp parallel &
            !$omp& private(j,k,tj,tk,head,group,cr,er,ar,c,b,q,rj,rk,rjk,rjk2, &
            !$omp& cosine,uj,uk,ujk,dcj,dck,qj,qk,qjk,dqj,dqk,dqjk,a,da,ca,cv,dj,dk,fc,dfc,ex,radial,rgradj,rgradk, &
            !$omp& last_group,last_cr,last_er,last_ar)
            do j = offsets(row),offsets(row+1)-1
                rj = geom(4,j)
                if (rj <= eps .or. rj > rcmax) cycle
                tj = local_species(species(indices(j)),s)
                if (tj < 1) cycle
                uj = geom(1:3,j)
                do k = j+1,offsets(row+1)-1
                    rk = geom(4,k)
                    if (rk <= eps .or. rk > rcmax) cycle
                    tk = local_species(species(indices(k)),s)
                    if (tk < 1) cycle
                    head = g4_first(tj,tk,row)
                    if (head == 0) cycle
                    ujk = dr(:,k)-dr(:,j); rjk2 = sum(ujk**2)
                    if (rjk2 <= eps**2 .or. rjk2 >= rcmax**2) cycle
                    rjk = sqrt(rjk2); ujk = ujk/rjk
                    uk = geom(1:3,k)
                    cosine = max(-1.0_real64,min(1.0_real64,sum(uj*uk)))
                    dcj = (uk-cosine*uj)/rj; dck = (uj-cosine*uk)/rk
                    ! Only the matching species-pair list is visited. Scalar
                    ! caches reuse consecutive radial/angular groups without
                    ! writing per-pair intermediates into global workspace.
                    last_group = 0; last_cr = 0; last_er = 0; last_ar = 0
                    !$omp do schedule(static)
                    do q = features(20,head,s),features(20,head,s)+features(21,head,s)-1
                        b = features(19,q,s)
                        if (rj <= fp(1,b,s) .and. rk <= fp(1,b,s) .and. rjk <= fp(1,b,s)) then
                            group = features(7,b,s)
                            if (group /= last_group) then
                                cr = features(15,b,s); er = features(16,b,s)
                                if (cr /= last_cr) then
                                    fc = target_cutoff_value(rjk,fp(1,b,s),features(4,b,s),fp(7,b,s))
                                    dfc = target_cutoff_derivative(rjk,fp(1,b,s),features(4,b,s),fp(7,b,s))
                                    last_cr = cr
                                end if
                                if (er /= last_er) then
                                    ex = exp(-fp(2,b,s)*(rjk-fp(3,b,s))**2)
                                    last_er = er
                                end if
                                qj = radial_cache(j,group,1); qk = radial_cache(k,group,1)
                                dqj = radial_cache(j,group,2); dqk = radial_cache(k,group,2)
                                qjk = fc*ex
                                dqjk = (dfc-2*fp(2,b,s)*(rjk-fp(3,b,s))*fc)*ex
                                radial = qj*qk*qjk
                                rgradj = dqj*uj*qk*qjk-qj*qk*dqjk*ujk
                                rgradk = qj*dqk*uk*qjk+qj*qk*dqjk*ujk
                                last_group = group
                            end if
                            ar = features(17,b,s)
                            if (ar /= last_ar) then
                                call angular_power(cosine,fp(4,b,s),fp(5,b,s),features(5,b,s), &
                                    0.5_real64*fp(5,b,s)*fp(4,b,s),a,da)
                                last_ar = ar
                            end if
                            g(row,b) = g(row,b)+2*a*radial
                            ca = 2*da*radial; cv = 2*a
                            dj = ca*dcj+cv*rgradj
                            dk = ca*dck+cv*rgradk
                            do c = 1,3
                                jacobian(c,b,j) = jacobian(c,b,j)+dj(c)
                                jacobian(c,b,k) = jacobian(c,b,k)+dk(c)
                            end do
                        end if
                    end do
                    !$omp end do nowait
                end do
            end do
            !$omp end parallel
        end do
        !$omp end target teams distribute
    end subroutine

    subroutine generic_values(device,nrw,meta,nodes,features,fp,local_species,species,centers,offsets,indices,geom,radial_cache,g)
        integer, intent(in) :: device,nrw
        integer, contiguous, intent(in) :: meta(:,:),nodes(:,:),features(:,:,:),local_species(:,:)
        integer, contiguous, intent(in) :: species(:),centers(:),offsets(:),indices(:)
        real(real64), contiguous, intent(in) :: fp(:,:,:),radial_cache(:,:,:),geom(:,:)
        real(real64), contiguous, intent(inout) :: g(:,:)
        integer :: row,b,s,j,k,tj,tk,kind,t1,t2,group,bb
        real(real64) :: total,v,dv,gradient(3),total12,v12,dv12
        !$omp target teams distribute parallel do collapse(2) device(device) if(device /= omp_get_initial_device()) &
        !$omp& map(alloc:meta,nodes,features,fp,local_species,species,centers,offsets,indices,geom,radial_cache,g) &
        !$omp& private(s,j,k,tj,tk,kind,t1,t2,group,bb,total,v,dv,gradient,total12,v12,dv12)
        do row = 1,nrw
            do b = 1,size(features,2)
                s = species(centers(row))
                if (meta(10,s) /= 1 .or. b > nodes(1,s)) cycle
                kind = features(1,b,s); t1 = features(2,b,s); t2 = features(3,b,s)
                if (kind == 4 .or. kind == 7) cycle
                if (features(10,b,s) /= 0 .and. features(10,b,s) /= b) cycle
                group = features(7,b,s)
                total = 0; total12 = 0
                do j = offsets(row),offsets(row+1)-1
                    tj = local_species(species(indices(j)),s)
                    if (kind /= 4 .and. kind /= 5) then
                        if (tj /= t1) cycle
                        if (kind == 6) then
                            call generic_lj(geom(4,j),features(4,b,s),fp(:,b,s),.false.,v,v12,dv,dv12)
                            total12 = total12+v12
                        else
                            v = generic_radial_value(kind,geom(4,j),features(4,b,s),fp(:,b,s))
                        end if
                        total = total+v
                    else
                        if (tj /= t1 .and. tj /= t2) cycle
                        do k = j+1,offsets(row+1)-1
                            tk = local_species(species(indices(k)),s)
                            if (.not. ((tj == t1 .and. tk == t2) .or. (tj == t2 .and. tk == t1))) cycle
                            v = generic_pair_value(features(:,b,s),fp(:,b,s),geom(1:3,j),geom(1:3,k), &
                                geom(4,j),geom(4,k),radial_cache(j,group,1),radial_cache(k,group,1))
                            total = total+v
                        end do
                    end if
                end do
                g(row,b) = total
                ! This path has one distinct angular power per group. Identical
                ! inputs still have separate NN weights, but share their value.
                bb = features(11,b,s)
                do while (bb /= 0)
                    g(row,bb) = total
                    bb = features(11,bb,s)
                end do
                if (kind == 6) g(row,b+1) = total12
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine

    ! Share pair geometry/radials within each angular group. Integer power sums
    ! remain local to the team, then each descriptor is written once. Stripping
    ! OpenMP directives gives the identical serial numerical loop.
    subroutine generic_values_grouped(device,nrw,meta,nodes,features,fp,local_species,species,centers,offsets, &
            indices,geom,radial_cache,g)
        integer, intent(in) :: device,nrw
        integer, contiguous, intent(in) :: meta(:,:),nodes(:,:),features(:,:,:),local_species(:,:)
        integer, contiguous, intent(in) :: species(:),centers(:),offsets(:),indices(:)
        real(real64), contiguous, intent(in) :: fp(:,:,:),radial_cache(:,:,:),geom(:,:)
        real(real64), contiguous, intent(inout) :: g(:,:)
        integer :: row,b,s,j,k,tj,tk,kind,t1,t2,group,bb,degree,d
        real(real64) :: total,v,dv,gradient(3),qjk,dqjk,ujk(3),cosine,total12,v12,dv12,t,power,hvalues(0:16)
        !$omp target teams distribute collapse(2) thread_limit(32) device(device) if(device /= omp_get_initial_device()) &
        !$omp& map(alloc:meta,nodes,features,fp,local_species,species,centers,offsets,indices,geom,radial_cache,g) &
        !$omp& private(s,j,k,tj,tk,kind,t1,t2,group,bb,degree,d,total,v,dv,gradient,qjk,dqjk,ujk,cosine, &
        !$omp& total12,v12,dv12,t,power,hvalues)
        do row = 1,nrw
            do b = 1,size(features,2)
                s = species(centers(row))
                if (meta(10,s) /= 1 .or. b > nodes(1,s)) cycle
                kind = features(1,b,s); t1 = features(2,b,s); t2 = features(3,b,s)
                if (kind == 4 .or. kind == 7) cycle
                if (features(10,b,s) /= 0 .and. features(10,b,s) /= b) cycle
                degree = features(12,b,s)
                hvalues = 0
                group = features(7,b,s)
                total = 0; total12 = 0
                !$omp parallel do private(j,k,tj,tk,d,v,dv,v12,dv12,qjk,dqjk,ujk,cosine,t,power) &
                !$omp& reduction(+:total,total12,hvalues)
                do j = offsets(row),offsets(row+1)-1
                    tj = local_species(species(indices(j)),s)
                    if (kind /= 4 .and. kind /= 5) then
                        if (tj /= t1) cycle
                        if (kind == 6) then
                            call generic_lj(geom(4,j),features(4,b,s),fp(:,b,s),.false.,v,v12,dv,dv12)
                            total12 = total12+v12
                        else
                            v = generic_radial_value(kind,geom(4,j),features(4,b,s),fp(:,b,s))
                        end if
                        total = total+v
                    else
                        if (tj /= t1 .and. tj /= t2) cycle
                        do k = j+1,offsets(row+1)-1
                            tk = local_species(species(indices(k)),s)
                            if (.not. ((tj == t1 .and. tk == t2) .or. (tj == t2 .and. tk == t1))) cycle
                            call generic_pair_geometry(features(:,b,s),fp(:,b,s),geom(1:3,j),geom(1:3,k), &
                                geom(4,j),geom(4,k),cosine,qjk,dqjk,ujk,.false.)
                            if (qjk == 0) cycle
                            v = 2*radial_cache(j,group,1)*radial_cache(k,group,1)*qjk
                            if (degree > 0) then
                                t = 0.5_real64*(1+fp(4,b,s)*cosine); power = 1
                                do d=1,degree
                                    power = power*t
                                    hvalues(d) = hvalues(d)+v*power
                                end do
                            else
                                total = total+v*angular_value(cosine,fp(4,b,s),fp(5,b,s),features(5,b,s))
                            end if
                        end do
                    end if
                end do
                !$omp end parallel do
                if (features(10,b,s) == b) then
                    bb = b
                    do while (bb /= 0)
                        if (degree > 0) then
                            g(row,bb) = hvalues(features(5,bb,s))
                        else
                            g(row,bb) = total
                        end if
                        bb = features(11,bb,s)
                    end do
                else
                    g(row,b) = total
                end if
                if (kind == 6) g(row,b+1) = total12
            end do
        end do
        !$omp end target teams distribute
    end subroutine

    subroutine generic_forces(device,nedges,meta,nodes,features,fp,local_species,species,centers,offsets, &
                              indices,geom,edge_row,radial_cache,jacobian,g,edge_force)
        integer, intent(in) :: device,nedges
        integer, contiguous, intent(in) :: meta(:,:),nodes(:,:),features(:,:,:),local_species(:,:)
        integer, contiguous, intent(in) :: species(:),centers(:),offsets(:),indices(:),edge_row(:)
        real(real64), contiguous, intent(in) :: fp(:,:,:),radial_cache(:,:,:),jacobian(:,:,:),geom(:,:),g(:,:)
        real(real64), contiguous, intent(inout) :: edge_force(:,:)
        integer :: row,b,s,j,k,tj,tk,kind,t1,t2,group,degree,d,c,bb
        real(real64) :: f(3),v,dv,gradient(3),gradient_k(3),v12,dv12,radial_force,hcoeff(0:16)
        logical :: pair_once
        ! CPU serial builds benefit from evaluating both sides once. GPU edge
        ! ownership avoids contended force atomics; all scalar formulas are shared.
        pair_once = device == omp_get_initial_device() .and. any(features(1,:,:) == 5)
        if (pair_once) then
            do j=1,nedges
                row=edge_row(j); s=species(centers(row))
                if (meta(10,s) == 1) edge_force(:,j)=0
            end do
        end if
        !$omp target teams distribute parallel do device(device) if(device /= omp_get_initial_device()) &
        !$omp& firstprivate(pair_once) map(alloc:meta,nodes,features,fp,local_species,species,centers,offsets,indices, &
        !$omp& geom,edge_row,radial_cache,jacobian,g,edge_force) &
        !$omp& private(row,b,s,k,tj,tk,kind,t1,t2,group,degree,d,c,bb,f,v,dv,gradient,gradient_k,v12,dv12,radial_force,hcoeff)
        do j = 1,nedges
            row = edge_row(j); s = species(centers(row))
            if (meta(10,s) /= 1) cycle
            tj = local_species(species(indices(j)),s); f = 0; radial_force = 0
            do b = 1,nodes(1,s)
                kind = features(1,b,s); t1 = features(2,b,s); t2 = features(3,b,s)
                if (kind == 7) cycle
                if (kind == 4) then
                    f = f+g(row,b)*jacobian(:,b,j)
                    cycle
                end if
                if (features(10,b,s) /= 0 .and. features(10,b,s) /= b) cycle
                group = features(7,b,s)
                if (kind /= 4 .and. kind /= 5) then
                    if (tj /= t1) cycle
                    if (kind == 6) then
                        call generic_lj(geom(4,j),features(4,b,s),fp(:,b,s),.true.,v,v12,dv,dv12)
                        radial_force = radial_force+g(row,b)*dv+g(row,b+1)*dv12
                    else
                        call generic_radial(kind,geom(4,j),features(4,b,s),fp(:,b,s),v,dv)
                        radial_force = radial_force+g(row,b)*dv
                    end if
                else
                    if (tj /= t1 .and. tj /= t2) cycle
                    degree = features(12,b,s)
                    hcoeff(:degree) = 0
                    if (degree > 0) then
                        bb=b
                        do while (bb /= 0)
                            if (features(13,bb,s) == bb) hcoeff(features(5,bb,s))=g(row,bb)
                            bb=features(11,bb,s)
                        end do
                    else
                        hcoeff(0)=g(row,b)
                    end if
                    do k = merge(j+1,offsets(row),pair_once),offsets(row+1)-1
                        if (k == j) cycle
                        tk = local_species(species(indices(k)),s)
                        if (.not. ((tj == t1 .and. tk == t2) .or. (tj == t2 .and. tk == t1))) cycle
                        call generic_pair_contracted(features(:,b,s),fp(:,b,s),geom(1:3,j),geom(1:3,k), &
                            geom(4,j),geom(4,k),radial_cache(j,group,1),radial_cache(k,group,1), &
                            radial_cache(j,group,2),radial_cache(k,group,2),hcoeff,degree,pair_once,gradient,gradient_k)
                        f = f+gradient
                        if (pair_once) then
                            do c=1,3
                                !$omp atomic update
                                edge_force(c,k) = edge_force(c,k)+gradient_k(c)
                            end do
                        end if
                    end do
                end if
            end do
            f = f+radial_force*geom(1:3,j)
            if (pair_once) then
                do c=1,3
                    !$omp atomic update
                    edge_force(c,j) = edge_force(c,j)+f(c)
                end do
            else
                edge_force(:,j) = f
            end if
        end do
        !$omp end target teams distribute parallel do
    end subroutine
end module
