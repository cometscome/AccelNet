    ! Decompose descriptor components only; keep the original NN intact.
    ! Each component uses the ordinary common descriptor/adjoint kernels.
    recursive subroutine initialize_composite(self,model,device,mode,status,message,use_host,g5_mode)
        class(target_model), intent(inout) :: self
        class(predictor_model), intent(in) :: model
        integer, optional, intent(in) :: device,mode,g5_mode
        integer, optional, intent(out) :: status
        character(len=*), optional, intent(out) :: message
        logical, optional, intent(in) :: use_host
        type(predictor_model) :: assembled,part
        type(descriptor_model) :: empty
        type(atomic_network) :: blank
        type(behler_config) :: radial
        integer :: ns,s,p,k,nc,nl,nb,np,d,j,st
        character(len=256) :: why
        ns=size(model%networks); assembled=model; np=0
        do s=1,ns
            nc=0; nl=0; nb=0
            if (allocated(model%setups(s)%model%chebyshev)) nc=size(model%setups(s)%model%chebyshev)
            if (allocated(model%setups(s)%model%lj)) nl=size(model%setups(s)%model%lj)
            if (allocated(model%setups(s)%model%behler)) nb=size(model%setups(s)%model%behler)
            np=max(np,nc+nl+nb)
            ! A shape-only descriptor shell validates/packs the full original NN.
            ! It is never evaluated: the component kernels populate its inputs.
            d=model%setups(s)%model%num_descriptors()
            assembled%setups(s)%model=empty
            call initialize_behler_config(radial,ns)
            do j=1,d
                call add_g1(radial,1,1.0_real64)
            end do
            call add_behler(assembled%setups(s)%model,radial)
        end do
        call self%initialize(assembled,device,mode,st,why,use_host,g5_mode)
        if (st/=0) then
            call target_failure(trim(why),status,message)
            return
        end if
        allocate(self%parts(np),self%part_offsets(ns,np),self%part_widths(ns,np))
        self%part_offsets=0; self%part_widths=0
        do p=1,np
            part=model
            do s=1,ns
                nc=0; nl=0; nb=0
                if (allocated(model%setups(s)%model%chebyshev)) nc=size(model%setups(s)%model%chebyshev)
                if (allocated(model%setups(s)%model%lj)) nl=size(model%setups(s)%model%lj)
                if (allocated(model%setups(s)%model%behler)) nb=size(model%setups(s)%model%behler)
                part%setups(s)%model=empty
                if (p<=nc) then
                    call add_chebyshev(part%setups(s)%model,model%setups(s)%model%chebyshev(p)%config)
                    self%part_offsets(s,p)=model%setups(s)%model%chebyshev(p)%output_offset
                else if (p<=nc+nl) then
                    k=p-nc
                    call add_lj(part%setups(s)%model,model%setups(s)%model%lj(k)%config)
                    self%part_offsets(s,p)=model%setups(s)%model%lj(k)%output_offset
                else if (p<=nc+nl+nb) then
                    k=p-nc-nl
                    call add_behler(part%setups(s)%model,model%setups(s)%model%behler(k)%config)
                    self%part_offsets(s,p)=model%setups(s)%model%behler(k)%output_offset
                else
                    call initialize_behler_config(radial,ns)
                    call add_g1(radial,1,1.0_real64)
                    call add_behler(part%setups(s)%model,radial)
                end if
                d=part%setups(s)%model%num_descriptors()
                if (p<=nc+nl+nb) self%part_widths(s,p)=d
                ! The child NN is never executed; only its descriptor dimensions
                ! and zero energy normalization are used for workspace packing.
                part%networks(s)=blank
                part%networks(s)%nlayers=2; part%networks(s)%maxnodes=d
                part%networks(s)%nodes=[d,1]; part%networks(s)%activation=[0]
                part%networks(s)%weight_offsets=[0]
                allocate(part%networks(s)%weights(d+1),part%networks(s)%descriptor_shift(d), &
                    part%networks(s)%descriptor_scale(d),part%networks(s)%atomic_references(ns))
                part%networks(s)%weights=0; part%networks(s)%descriptor_shift=0
                part%networks(s)%descriptor_scale=1; part%networks(s)%atomic_references=0
            end do
            call self%parts(p)%initialize(part,device,mode,st,why,use_host,g5_mode)
            if (st/=0) then
                call self%release()
                call target_failure(trim(why),status,message)
                return
            end if
        end do
        if (present(status)) status=0
        if (present(message)) message=''
    end subroutine

    subroutine execute_composite(model,work,nrw,natoms,nedges,stages,energy_only,descriptors_only)
        type(target_model), intent(in) :: model
        type(target_workspace), intent(inout) :: work
        integer, intent(in) :: nrw,natoms,nedges
        real(real64), intent(out) :: stages(3)
        logical, optional, intent(in) :: energy_only
        logical, intent(in) :: descriptors_only
        real(real64) :: part_stages(3)
        integer :: p,q
        logical :: only_energy
        only_energy=.false.
        if (present(energy_only)) only_energy=energy_only
        stages=0
        if (allocated(work%parts)) then
            if (size(work%parts)/=size(model%parts)) then
                do q=1,size(work%parts)
                    call work%parts(q)%release()
                end do
                deallocate(work%parts)
            end if
        end if
        if (.not.allocated(work%parts)) allocate(work%parts(size(model%parts)))
        call clear_composite(model%device,nrw,natoms,work%g,work%forces,work%virial)
        do p=1,size(model%parts)
            call prepare_workspace(model%parts(p),work%parts(p),nrw,natoms,nedges)
            call copy_composite_csr(model%device,nrw,natoms,nedges,work%species,work%centers,work%offsets, &
                work%indices,work%dr,work%parts(p)%species,work%parts(p)%centers,work%parts(p)%offsets, &
                work%parts(p)%indices,work%parts(p)%dr)
            call execute_phase(model%parts(p),work%parts(p),nrw,natoms,nedges,part_stages,phase=1)
            stages=stages+part_stages
            call transfer_component(model%device,nrw,p,model%part_offsets,model%part_widths,work%species, &
                work%centers,work%g,work%parts(p)%g,.false.)
        end do
        if (descriptors_only) return
        call execute_phase(model,work,nrw,natoms,nedges,part_stages,energy_only,phase=2)
        stages=stages+part_stages
        if (only_energy) return
        do p=1,size(model%parts)
            call transfer_component(model%device,nrw,p,model%part_offsets,model%part_widths,work%species, &
                work%centers,work%g,work%parts(p)%g,.true.)
            call execute_phase(model%parts(p),work%parts(p),nrw,natoms,nedges,part_stages,phase=3)
            stages=stages+part_stages
            call add_component_forces(model%device,natoms,work%forces,work%virial, &
                work%parts(p)%forces,work%parts(p)%virial)
        end do
    end subroutine

    subroutine clear_composite(device,nrw,natoms,g,forces,virial)
        integer, intent(in) :: device,nrw,natoms
        real(real64), contiguous, intent(inout) :: g(:,:),forces(:,:),virial(:,:)
        integer :: i
        !$omp target teams distribute parallel do device(device) if(target:device/=omp_get_initial_device()) &
        !$omp& map(alloc:g,forces,virial)
        do i=1,max(nrw,natoms)
            if (i<=nrw) g(i,:)=0
            if (i<=natoms) forces(:,i)=0
            if (i==1) virial=0
        end do
        !$omp end target teams distribute parallel do
    end subroutine

    subroutine copy_composite_csr(device,nrw,natoms,nedges,sp,ce,off,ind,dr,sp2,ce2,off2,ind2,dr2)
        integer, intent(in) :: device,nrw,natoms,nedges
        integer, contiguous, intent(in) :: sp(:),ce(:),off(:),ind(:)
        integer, contiguous, intent(inout) :: sp2(:),ce2(:),off2(:),ind2(:)
        real(real64), contiguous, intent(in) :: dr(:,:)
        real(real64), contiguous, intent(inout) :: dr2(:,:)
        integer :: i
        !$omp target teams distribute parallel do device(device) if(target:device/=omp_get_initial_device()) &
        !$omp& map(alloc:sp,ce,off,ind,dr,sp2,ce2,off2,ind2,dr2)
        do i=1,max(nrw+1,natoms,nedges)
            if (i<=natoms) sp2(i)=sp(i)
            if (i<=nrw) ce2(i)=ce(i)
            if (i<=nrw+1) off2(i)=off(i)
            if (i<=nedges) then
                ind2(i)=ind(i); dr2(:,i)=dr(:,i)
            end if
        end do
        !$omp end target teams distribute parallel do
    end subroutine

    subroutine transfer_component(device,nrw,p,offsets,widths,species,centers,g,part,adjoint)
        integer, intent(in) :: device,nrw,p
        integer, contiguous, intent(in) :: offsets(:,:),widths(:,:),species(:),centers(:)
        real(real64), contiguous, intent(inout) :: g(:,:),part(:,:)
        logical, intent(in) :: adjoint
        integer :: row,s,b
        !$omp target teams distribute parallel do device(device) if(target:device/=omp_get_initial_device()) &
        !$omp& map(to:offsets,widths) map(alloc:species,centers,g,part) private(s,b)
        do row=1,nrw
            s=species(centers(row))
            if (adjoint) part(row,:)=0
            do b=1,widths(s,p)
                if (adjoint) then
                    part(row,b)=g(row,offsets(s,p)+b)
                else
                    g(row,offsets(s,p)+b)=part(row,b)
                end if
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine

    subroutine add_component_forces(device,natoms,forces,virial,part_force,part_virial)
        integer, intent(in) :: device,natoms
        real(real64), contiguous, intent(inout) :: forces(:,:),virial(:,:)
        real(real64), contiguous, intent(in) :: part_force(:,:),part_virial(:,:)
        integer :: i
        !$omp target teams distribute parallel do device(device) if(target:device/=omp_get_initial_device()) &
        !$omp& map(alloc:forces,virial,part_force,part_virial)
        do i=1,natoms
            forces(:,i)=forces(:,i)+part_force(:,i)
            if (i==1) virial=virial+part_virial
        end do
        !$omp end target teams distribute parallel do
    end subroutine

    subroutine download_descriptors(device,nrw,g,values)
        integer, intent(in) :: device,nrw
        real(real64), contiguous, intent(inout) :: g(:,:)
        real(real64), intent(out) :: values(:,:)
        !$omp target update from(g) device(device) if(device/=omp_get_initial_device())
        values=g(:nrw,:size(values,2))
    end subroutine
