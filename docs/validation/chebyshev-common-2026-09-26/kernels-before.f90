! All arguments are plain contiguous arrays. No derived-type deep mapping,
! device allocation, or descriptor Jacobian is needed by the kernels.
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
                                energies, forces, virial, g, values, deriv, delta, moments, powers, use_moment, &
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
        real(real64), contiguous, intent(inout) :: moments(:,:,:), powers(:,:,:)
        integer, contiguous, intent(inout) :: use_moment(:), edge_row(:)
        real(real64), contiguous, intent(inout) :: geom(:,:), edge_force(:,:)
        real(real64), intent(out) :: stages(3)
        integer :: row, s, nr, na, dim, multi, version, ct, j, k, b, l, i, o, nin, nout, a, c, target
        integer :: entry, q, ax, ay, az, angular_neighbors, nm
        real(real64) :: rj, rk, fcj, fck, dfcj, sj, sk, cosine, x, t0, t1, t2, dt0, dt1, dt2
        real(real64) :: v, dv, z, rc, ac, alpha, xscale, vc, dc, coeff, radial, fj(3), center_f(3), w(3,3)
        real(real64) :: self0, self1, mono, u(3), grad(3), dm(3), correction, started
        ! All arrays have persistent mappings owned by target_workspace.
        ! Structured references here neither copy data nor allocate device memory.
        !$omp target data device(device) if(device /= omp_get_initial_device()) &
        !$omp& map(alloc: features,feature_params,local_species) &
        !$omp& map(alloc: meta, nodes, acts, woffset, weights, params, shift, scale, spin, mp, multiplicity, polynomial) &
        !$omp& map(alloc: species, centers, offsets, indices, dr, energies, forces, virial, &
        !$omp& g, values, deriv, delta, moments, powers, use_moment, geom, edge_row, edge_force)
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
        if (any(meta(10,:) == 1)) call generic_values(device,nrw,meta,nodes,features,feature_params, &
            local_species,species,centers,offsets,indices,geom,g)
        stages(1) = omp_get_wtime()-started
        started = omp_get_wtime()
        !$omp target teams distribute parallel do device(device) if(device /= omp_get_initial_device()) thread_limit(32) &
        !$omp& private(s,dim,l,nin,nout,o,j,i,z,v,b,q,nr,na,multi,entry,self0,self1)
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
        !$omp& entry,q,ax,ay,az,nm,mono,u,grad,dm,correction,v)
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
                    fj = radial*dr(:,j)/rj
                end if
                if (rj <= ac .and. rj >= eps) then
                    fcj = geom(5,j); dfcj = geom(6,j)
                    if (use_moment(row) == 1) then
                        u = geom(1:3,j)
                        v = 0; dm = 0
                        do entry = 1, nm
                            ax = mp(1,entry,s); ay = mp(2,entry,s); az = mp(3,entry,s)
                            mono = powers(j,ax+1,1)*powers(j,ay+1,2)*powers(j,az+1,3)
                            grad = 0
                            if (ax > 0) grad(1) = ax*powers(j,ax,1)*powers(j,ay+1,2)*powers(j,az+1,3)
                            if (ay > 0) grad(2) = ay*powers(j,ax+1,1)*powers(j,ay,2)*powers(j,az+1,3)
                            if (az > 0) grad(3) = az*powers(j,ax+1,1)*powers(j,ay+1,2)*powers(j,az,3)
                            coeff = moments(row,entry,1)+sj*moments(row,entry,2)
                            v = v+coeff*mono
                            dm = dm+coeff*grad
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
        if (any(meta(10,:) == 1)) call generic_forces(device,nedges,meta,nodes,features,feature_params, &
            local_species,species,centers,offsets,indices,geom,edge_row,g,edge_force)
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
    subroutine generic_values(device,nrw,meta,nodes,features,fp,local_species,species,centers,offsets,indices,geom,g)
        integer, intent(in) :: device,nrw
        integer, contiguous, intent(in) :: meta(:,:),nodes(:,:),features(:,:,:),local_species(:,:)
        integer, contiguous, intent(in) :: species(:),centers(:),offsets(:),indices(:)
        real(real64), contiguous, intent(in) :: fp(:,:,:),geom(:,:)
        real(real64), contiguous, intent(inout) :: g(:,:)
        integer :: row,b,s,j,k,tj,tk,kind,t1,t2
        real(real64) :: total,v,dv,gradient(3)
        !$omp target teams distribute parallel do collapse(2) device(device) if(device /= omp_get_initial_device()) &
        !$omp& map(alloc:meta,nodes,features,fp,local_species,species,centers,offsets,indices,geom,g) &
        !$omp& private(s,j,k,tj,tk,kind,t1,t2,total,v,dv,gradient)
        do row = 1,nrw
            do b = 1,size(features,2)
                s = species(centers(row))
                if (meta(10,s) /= 1 .or. b > nodes(1,s)) cycle
                kind = features(1,b,s); t1 = features(2,b,s); t2 = features(3,b,s)
                total = 0
                do j = offsets(row),offsets(row+1)-1
                    tj = local_species(species(indices(j)),s)
                    if (kind /= 4 .and. kind /= 5) then
                        if (tj /= t1) cycle
                        call generic_radial(kind,geom(4,j),features(4,b,s),fp(:,b,s),v,dv)
                        total = total+v
                    else
                        if (tj /= t1 .and. tj /= t2) cycle
                        do k = j+1,offsets(row+1)-1
                            tk = local_species(species(indices(k)),s)
                            if (.not. ((tj == t1 .and. tk == t2) .or. (tj == t2 .and. tk == t1))) cycle
                            call generic_pair(features(:,b,s),fp(:,b,s),geom(1:3,j),geom(1:3,k), &
                                geom(4,j),geom(4,k),v,gradient)
                            total = total+v
                        end do
                    end if
                end do
                g(row,b) = total
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine

    subroutine generic_forces(device,nedges,meta,nodes,features,fp,local_species,species,centers,offsets, &
                              indices,geom,edge_row,g,edge_force)
        integer, intent(in) :: device,nedges
        integer, contiguous, intent(in) :: meta(:,:),nodes(:,:),features(:,:,:),local_species(:,:)
        integer, contiguous, intent(in) :: species(:),centers(:),offsets(:),indices(:),edge_row(:)
        real(real64), contiguous, intent(in) :: fp(:,:,:),geom(:,:),g(:,:)
        real(real64), contiguous, intent(inout) :: edge_force(:,:)
        integer :: row,b,s,j,k,tj,tk,kind,t1,t2
        real(real64) :: f(3),v,dv,gradient(3)
        !$omp target teams distribute parallel do device(device) if(device /= omp_get_initial_device()) &
        !$omp& map(alloc:meta,nodes,features,fp,local_species,species,centers,offsets,indices,geom,edge_row,g,edge_force) &
        !$omp& private(row,b,s,k,tj,tk,kind,t1,t2,f,v,dv,gradient)
        do j = 1,nedges
            row = edge_row(j); s = species(centers(row))
            if (meta(10,s) /= 1) cycle
            tj = local_species(species(indices(j)),s); f = 0
            do b = 1,nodes(1,s)
                kind = features(1,b,s); t1 = features(2,b,s); t2 = features(3,b,s)
                if (kind /= 4 .and. kind /= 5) then
                    if (tj /= t1) cycle
                    call generic_radial(kind,geom(4,j),features(4,b,s),fp(:,b,s),v,dv)
                    f = f+g(row,b)*dv*geom(1:3,j)
                else
                    if (tj /= t1 .and. tj /= t2) cycle
                    do k = offsets(row),offsets(row+1)-1
                        if (k == j) cycle
                        tk = local_species(species(indices(k)),s)
                        if (.not. ((tj == t1 .and. tk == t2) .or. (tj == t2 .and. tk == t1))) cycle
                        call generic_pair(features(:,b,s),fp(:,b,s),geom(1:3,j),geom(1:3,k), &
                            geom(4,j),geom(4,k),v,gradient)
                        f = f+g(row,b)*gradient
                    end do
                end if
            end do
            edge_force(:,j) = f
        end do
        !$omp end target teams distribute parallel do
    end subroutine
end module
