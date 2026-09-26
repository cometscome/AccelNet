from pathlib import Path
root=Path('/tmp/accelnet-gpu-stage4/mapping-check'); root.mkdir(exist_ok=True)
p=Path('/home/nagai/AccelNetGPU/AccelNet-clone/AccelNetPredictor/src/accelnet_batch_target.f90')
s=p.read_text().replace('    use iso_fortran_env, only: real64','    use iso_fortran_env, only: real64\n    use iso_c_binding, only: c_ptr,c_loc\n    use omp_lib, only: omp_target_is_present')
for name,names in [('map_model_arrays','meta,nodes,acts,woffset,weights,params,shift,scale,spin,mp,multiplicity,polynomial'),('map_buffers','g,values,deriv,delta,moments,powers,species,centers,offsets,indices,use_moment,dr,energies,forces,virial,geom,edge_row,edge_force')]:
    start=s.index('    subroutine '+name+'('); end=s.index('    end subroutine',start)
    block=s[start:end].replace('contiguous, intent','contiguous, target, intent')
    calls='\n'.join('            call require_unmapped(c_loc('+n+'),device)' for n in names.split(','))
    block=block.replace('        end if',calls+'\n        end if')
    s=s[:start]+block+s[end:]
s=s.replace('end module','''    subroutine require_unmapped(ptr,device)
        type(c_ptr), value :: ptr
        integer :: device
        if (omp_target_is_present(ptr,device) /= 0) error stop 'mapping remained after delete'
    end subroutine
end module''')
(root/'accelnet_batch_target.f90').write_text(s)
