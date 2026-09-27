! Common edge force arithmetic, included by center-owned and edge-owned schedules.
            row = edge_row(j); s = species(centers(row))
            if (meta(10,s) /= 1) cycle
            if (meta(11,s) > 0) then
                if (g5_active(row,2*size(features,2)+1) == 1) then
                    edge_force(:,j) = 0
                    cycle
                end if
            end if
            tj = local_species(species(indices(j)),s); f = 0; radial_force = 0
            do b = 1,nodes(1,s)
                kind = features(1,b,s); t1 = features(2,b,s); t2 = features(3,b,s)
                if (kind == 7) cycle
                if (meta(11,s) > 0) then
                    if (g5_active(row,b) == 1) cycle
                end if
                if (kind == 4 .or. kind >= 12) then
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
                    else if (kind == 2) then
                        radial_force = radial_force+g(row,b)*radial_cache(j,group,2)
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
