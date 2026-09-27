! Packed non-Chebyshev features. One column per output in the CPU model.
! Integer fields: kind (G1..G5=1..5, LJ6=6, LJ12=7), species pair,
! cutoff type, integer zeta, local species count, angular radial group, group representative.
! Fields 9..12: angular group, representative, next member, maximum polynomial degree.
! Fields 13/14: equal-power representative and next equal-power member.
! Real fields: Rc, eta, Rs, lambda, zeta, kappa, cutoff alpha.
module accelnet_target_descriptors
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet_descriptor_models, only: descriptor_model
    implicit none
    private
    public :: pack_generic_descriptors, PACKED_FEATURE_FIELDS, PACKED_PARAMETER_FIELDS
    integer, parameter :: PACKED_FEATURE_FIELDS = 31
    integer, parameter :: PACKED_PARAMETER_FIELDS = 19
contains
    subroutine pack_generic_descriptors(model, fi, fr, status, message, g5_mode)
        type(descriptor_model), intent(in) :: model
        integer, allocatable, intent(out) :: fi(:,:)
        real(real64), allocatable, intent(out) :: fr(:,:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer, optional, intent(in) :: g5_mode
        integer :: c,j,b,o,t,n,ct,ns,groups,k,r,last,q,zeta,ne,first,last_ext,nangle
        integer, allocatable :: ext_order(:)
        real(real64) :: binomial
        real(real64) :: alpha
        status = 1
        message = 'OpenMP target: invalid LJ/Behler descriptor metadata'
        n = model%num_outputs
        allocate(fi(PACKED_FEATURE_FIELDS,n),fr(PACKED_PARAMETER_FIELDS,n))
        fi = 0; fr = 0
        if (n < 1) return
        if (allocated(model%lj)) then
            do c = 1,size(model%lj)
                associate(cfg => model%lj(c)%config)
                o = model%lj(c)%output_offset; ns = cfg%num_species
                if (o < 0 .or. ns < 1 .or. o+2*ns > n) return
                do t = 1,ns
                    do j = 1,2
                        b = o+2*(t-1)+j
                        if (fi(1,b) /= 0) return
                        fi(:6,b) = [5+j,t,0,cfg%cutoff_type,0,ns]
                        fr(1,b) = cfg%radial_rc; fr(7,b) = cfg%cutoff_alpha
                    end do
                end do
                end associate
            end do
        end if
        if (allocated(model%behler)) then
            do c = 1,size(model%behler)
                associate(cfg => model%behler(c)%config)
                o = model%behler(c)%output_offset; ns = cfg%num_species
                ct = cfg%cutoff_type; alpha = cfg%cutoff_alpha
                if (allocated(cfg%g1)) then
                    do j = 1,size(cfg%g1)
                        associate(p => cfg%g1(j))
                        b = o+p%output
                        if (b < 1 .or. b > n) return
                        if (fi(1,b) /= 0) return
                        fi(:6,b) = [1,p%species,0,ct,0,ns]
                        fr(1,b) = p%rc; fr(7,b) = alpha
                        end associate
                    end do
                end if
                if (allocated(cfg%g2)) then
                    do j = 1,size(cfg%g2)
                        associate(p => cfg%g2(j))
                        b = o+p%output
                        if (b < 1 .or. b > n) return
                        if (fi(1,b) /= 0) return
                        fi(:6,b) = [2,p%species,0,ct,0,ns]
                        fr(1,b) = p%rc; fr(7,b) = alpha
                        fr(2,b) = p%eta; fr(3,b) = p%rs
                        end associate
                    end do
                end if
                if (allocated(cfg%g3)) then
                    do j = 1,size(cfg%g3)
                        associate(p => cfg%g3(j))
                        b = o+p%output
                        if (b < 1 .or. b > n) return
                        if (fi(1,b) /= 0) return
                        fi(:6,b) = [3,p%species,0,ct,0,ns]
                        fr(1,b) = p%rc; fr(7,b) = alpha
                        fr(6,b) = p%kappa
                        end associate
                    end do
                end if
                if (allocated(cfg%g4)) then
                    do j = 1,size(cfg%g4)
                        associate(p => cfg%g4(j))
                        b = o+p%output
                        if (b < 1 .or. b > n) return
                        if (fi(1,b) /= 0) return
                        fi(:6,b) = [4,p%species1,p%species2,ct,p%integer_zeta,ns]
                        fr(1,b) = p%rc; fr(7,b) = alpha
                        fr(2,b) = p%eta; fr(3,b) = p%rs
                        fr(4,b) = p%lambda; fr(5,b) = p%zeta
                        end associate
                    end do
                end if
                if (allocated(cfg%g5)) then
                    do j = 1,size(cfg%g5)
                        associate(p => cfg%g5(j))
                        b = o+p%output
                        if (b < 1 .or. b > n) return
                        if (fi(1,b) /= 0) return
                        fi(:6,b) = [5,p%species1,p%species2,ct,p%integer_zeta,ns]
                        fi(22,b) = cfg%g5_evaluation_mode
                        if (present(g5_mode)) fi(22,b) = g5_mode
                        fr(19,b) = cfg%maximum_angular_cutoff
                        if (fi(22,b) < 0 .or. fi(22,b) > 3) return
                        fr(1,b) = p%rc; fr(7,b) = alpha
                        fr(2,b) = p%eta; fr(3,b) = p%rs
                        fr(4,b) = p%lambda; fr(5,b) = p%zeta
                        end associate
                    end do
                end if
                if (allocated(cfg%extended)) then
                    do j=1,size(cfg%extended)
                        associate(p => cfg%extended(j))
                        b=o+p%output
                        if (b < 1 .or. b > n) return
                        if (fi(1,b) /= 0) return
                        fi(:6,b)=[p%kind,p%species1,p%species2,ct,0,ns]
                        fr(1:7,b)=p%p
                        end associate
                    end do
                end if
                end associate
            end do
        end if
        if (.not. all(ieee_is_finite(fr))) return
        do b = 1,n
            if (fi(1,b) >= 12) then
                if (.not. (fi(1,b) == 12 .or. fi(1,b) == 13 .or. &
                    (fi(1,b) >= 20 .and. fi(1,b) <= 25))) return
                if (fr(1,b) <= 0) return
                cycle ! add_extended already validates compact and weighted parameters.
            end if
            if (fi(1,b) < 1 .or. fi(1,b) > 7) return
            if (fi(2,b) < 1 .or. fi(2,b) > fi(6,b)) return
            if (fr(1,b) <= 0) return
            ct = fi(4,b); alpha = fr(7,b)
            if (ct < 0 .or. ct > 9 .or. alpha < 0 .or. alpha >= 1) return
            if (ct == 9 .and. alpha <= 0) return
            if (fi(1,b) == 4 .or. fi(1,b) == 5) then
                if (fi(3,b) < 1 .or. fi(3,b) > fi(6,b)) return
                if (abs(fr(4,b)) > 1 .or. fr(5,b) < 1) return
            end if
        end do
        ! G2, G4 and G5 have the same q(r). Cache its value/derivative once
        ! across species pairs and angular powers, including radial-only G2.
        groups = 0
        do b = 1,n
            if (fi(1,b) /= 2 .and. fi(1,b) /= 4 .and. fi(1,b) /= 5) cycle
            do k = 1,b-1
                if (fi(7,k) == 0) cycle
                if (fi(4,b) /= fi(4,k)) cycle
                if (any(fr(1:3,b) /= fr(1:3,k)) .or. fr(7,b) /= fr(7,k)) cycle
                fi(7:8,b) = fi(7:8,k)
                exit
            end do
            if (fi(7,b) /= 0) cycle
            groups = groups+1
            fi(7:8,b) = [groups,b]
        end do
        ! A G2-only radial group may need just one neighbor species. Field 3
        ! of its representative records that filter (zero means unrestricted).
        do b=1,n
            if (fi(1,b) /= 2 .or. fi(8,b) /= b) cycle
            fi(3,b)=fi(2,b)
            do k=1,n
                if (fi(7,k) /= fi(7,b)) cycle
                if (fi(1,k) /= 2 .or. fi(2,k) /= fi(2,b)) then
                    fi(3,b)=0
                    exit
                end if
            end do
        end do
        ! Extended descriptors share radial values across angular windows/species.
        do b=1,n
            if (fi(1,b) < 12) cycle
            do k=1,b-1
                if (fi(1,k) < 12) cycle
                if ((fi(1,b) >= 20) .neqv. (fi(1,k) >= 20)) cycle
                if (fi(4,b) /= fi(4,k)) cycle
                if (fi(1,b) >= 20) then
                    if (any(fr(1:2,b) /= fr(1:2,k)) .or. fr(5,b) /= fr(5,k)) cycle
                else
                    if (any(fr(1:3,b) /= fr(1:3,k)) .or. fr(7,b) /= fr(7,k)) cycle
                end if
                fi(7:8,b)=fi(7:8,k); exit
            end do
            if (fi(7,b) /= 0) cycle
            groups=groups+1; fi(7:8,b)=[groups,b]
        end do
        ! G5 integer moments: compact only radial groups actually requested.
        ! Exact integer powers 1..10 use a polynomial; fractional/near-integer
        ! and higher powers keep the direct formula. Auto retains the original
        ! per-component angular-neighbor threshold. Fields 22/23/24: mode/group/rep.
        groups = 0
        do b = 1,n
            if (fi(1,b) /= 5 .or. fi(22,b) == 1) cycle
            zeta = fi(5,b)
            if (zeta < 1 .or. zeta > 10) cycle
            if (fr(5,b) /= real(zeta,real64)) cycle
            do k = 1,b-1
                if (fi(23,k) == 0 .or. fi(7,k) /= fi(7,b)) cycle
                fi(23:24,b) = fi(23:24,k)
                exit
            end do
            if (fi(23,b) == 0) then
                groups = groups+1
                fi(23:24,b) = [groups,b]
                fi(25,groups) = b
            end if
            binomial = 1
            do q = 0,zeta
                fr(8+q,b) = 2.0_real64**(1-zeta)*binomial*fr(4,b)**q
                binomial = binomial*real(zeta-q,real64)/real(q+1,real64)
            end do
        end do
        ! Contract only features with identical non-angular factors and lambda.
        ! Positive integer zeta <= 16 shares a polynomial in t=(1+lambda*c)/2.
        ! Fractional/high zeta retains its exact power and only identical powers combine.
        groups = 0
        do b = 1,n
            if (fi(1,b) /= 4 .and. fi(1,b) /= 5) cycle
            r = 0
            do k = 1,b-1
                if (fi(10,k) /= k) cycle
                if (fi(1,b) /= fi(1,k) .or. fi(7,b) /= fi(7,k)) cycle
                if (fi(22,b) /= fi(22,k) .or. fr(19,b) /= fr(19,k)) cycle
                if ((fi(23,b) > 0) .neqv. (fi(23,k) > 0)) cycle
                if (minval(fi(2:3,b)) /= minval(fi(2:3,k))) cycle
                if (maxval(fi(2:3,b)) /= maxval(fi(2:3,k))) cycle
                if (fr(4,b) /= fr(4,k)) cycle
                if (fi(5,b) > 0 .and. fi(5,b) <= 16 .and. fr(5,b) == real(fi(5,b),real64) .and. fi(12,k) > 0) then
                    r = k
                else if (fi(12,k) == 0 .and. fr(5,b) == fr(5,k)) then
                    r = k
                end if
                if (r /= 0) exit
            end do
            if (r == 0) then
                groups = groups+1
                fi(9,b) = groups; fi(10,b) = b
                if (fi(5,b) > 0 .and. fi(5,b) <= 16 .and. fr(5,b) == real(fi(5,b),real64)) fi(12,b) = fi(5,b)
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
        ! A single distinct power needs no polynomial recurrence. Keep the
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
        ! Equal-power coefficients occupy their original first NN-gradient slot.
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
        ! G4 follows the established CPU value/Jacobian traversal. Reuse cutoff,
        ! exponential, angular and radial factors independently, then follow a
        ! species-pair list instead of scanning all descriptors for each pair.
        ! Fields 15..18: cutoff rep, exponential rep, angular rep (within species
        ! pair), and next G4 in species pair. Fields 19..21 pack indexed lists.
        do b = 1,n
            if (fi(1,b) /= 4) cycle
            fi(15:17,b) = b
            fr(8,b) = 0.5_real64*fr(5,b)*fr(4,b) ! invariant angular derivative prefactor
            do k = 1,b-1
                if (fi(1,k) /= 4) cycle
                if (fi(4,k) == fi(4,b) .and. fr(1,k) == fr(1,b) .and. fr(7,k) == fr(7,b)) &
                    fi(15,b) = fi(15,k)
                if (all(fr(2:3,k) == fr(2:3,b))) fi(16,b) = fi(16,k)
                if (minval(fi(2:3,k)) /= minval(fi(2:3,b)) .or. &
                    maxval(fi(2:3,k)) /= maxval(fi(2:3,b))) cycle
                if (all(fr(4:5,k) == fr(4:5,b))) fi(17,b) = fi(17,k)
                ! Ascending construction appends b to the previous pair member.
                if (fi(18,k) == 0) then
                    fi(18,k) = b
                end if
            end do
        end do
        ! Flatten each species-pair list. Static worksharing assigns each
        ! descriptor to the same thread for every pair of a center.
        groups = 0
        do b = 1,n
            if (fi(1,b) /= 4) cycle
            if (any(fi(18,:) == b)) cycle
            fi(20,b) = groups+1
            r = b
            do while (r /= 0)
                groups = groups+1
                fi(19,groups) = r
                r = fi(18,r)
            end do
            fi(21,b) = groups-fi(20,b)+1
        end do
        ! Extended angular ownership and scalar-cache ordering. Fields 26..31:
        ! sorted feature index, pair-list start/count, angular representative,
        ! weighted-list head (column 1), total angular count (column 1).
        ! Feature output indices and NN/scaling order are never permuted.
        if (any(fi(1,:) == 13 .or. fi(1,:) == 21 .or. fi(1,:) == 22 .or. fi(1,:) == 24 .or. fi(1,:) == 25)) then
            allocate(ext_order(n)); ne=0; nangle=0
            do b=1,n
                if (.not. extended_angular_kind(fi(1,b))) cycle
                if (fi(1,b) /= 13) then
                    fr(8,b)=2.0_real64/(fr(4,b)-fr(3,b))
                    fr(9,b)=0.5_real64*(fr(3,b)+fr(4,b))
                end if
                ne=ne+1; ext_order(ne)=b; fi(29,b)=0
                do k=1,b-1
                    if (.not. extended_angular_kind(fi(1,k))) cycle
                    if ((fi(1,b) == 13) .neqv. (fi(1,k) == 13)) cycle
                    if (fi(1,b) == 13) then
                        if (any(fr(4:5,b) /= fr(4:5,k))) cycle
                    else
                        if (any(fr(3:4,b) /= fr(3:4,k))) cycle
                        if (mod(nint(fr(5,b)),10) /= mod(nint(fr(5,k)),10)) cycle
                    end if
                    fi(29,b)=fi(29,k); exit
                end do
                if (fi(29,b) == 0) then
                    nangle=nangle+1; fi(29,b)=nangle
                end if
            end do
            ! Reuse the G4 cutoff-cache field for weighted angular members.
            do b=1,n
                if (fi(1,b) /= 13) cycle
                fi(15,b)=b
                do k=1,b-1
                    if (fi(1,k) /= 13 .or. fi(4,b) /= fi(4,k)) cycle
                    if (fr(1,b) /= fr(1,k) .or. fr(7,b) /= fr(7,k)) cycle
                    fi(15,b)=fi(15,k); exit
                end do
            end do
            do j=2,ne
                b=ext_order(j); k=j-1
                do while(k >= 1)
                    if (.not. extended_less(b,ext_order(k),fi)) exit
                    ext_order(k+1)=ext_order(k); k=k-1
                end do
                ext_order(k+1)=b
            end do
            fi(26,1:ne)=ext_order(1:ne); fi(31,1)=ne
            first=1
            do while(first <= ne)
                b=ext_order(first); last_ext=first
                do while(last_ext < ne)
                    k=ext_order(last_ext+1)
                    if (minval(fi(2:3,b)) /= minval(fi(2:3,k)) .or. &
                        maxval(fi(2:3,b)) /= maxval(fi(2:3,k))) exit
                    last_ext=last_ext+1
                end do
                fi(27:28,b)=[first,last_ext-first+1]
                if (fi(2,b) == 0) fi(30,1)=b
                ! Field 16 is the packed end of a constant-angular/narrow run.
                ! G4 uses that field for its separate exponential cache key.
                j=first
                do while(j <= last_ext)
                    b=ext_order(j); q=j
                    do while(q < last_ext)
                        k=ext_order(q+1)
                        if (fi(29,b) /= fi(29,k)) exit
                        if ((fi(1,b) == 13 .or. fi(1,b) == 21 .or. fi(1,b) == 24) .neqv. &
                            (fi(1,k) == 13 .or. fi(1,k) == 21 .or. fi(1,k) == 24)) exit
                        q=q+1
                    end do
                    fi(16,ext_order(j:q))=q
                    j=q+1
                end do
                first=last_ext+1
            end do
        end if
        status = 0; message = ''
    end subroutine
    pure logical function extended_angular_kind(kind) result(angular)
        integer, intent(in) :: kind
        angular=kind == 13 .or. kind == 21 .or. kind == 22 .or. kind == 24 .or. kind == 25
    end function

    pure logical function extended_less(b,k,fi) result(less)
        integer, intent(in) :: b,k,fi(:,:)
        integer :: left(5),right(5),i
        left=[minval(fi(2:3,b)),maxval(fi(2:3,b)),fi(29,b), &
            merge(1,0,fi(1,b) == 13 .or. fi(1,b) == 21 .or. fi(1,b) == 24),fi(7,b)]
        right=[minval(fi(2:3,k)),maxval(fi(2:3,k)),fi(29,k), &
            merge(1,0,fi(1,k) == 13 .or. fi(1,k) == 21 .or. fi(1,k) == 24),fi(7,k)]
        less=.false.
        do i=1,5
            if (left(i) == right(i)) cycle
            less=left(i)<right(i); return
        end do
    end function
end module
