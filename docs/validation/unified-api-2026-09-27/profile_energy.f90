program profile_energy
use iso_fortran_env, only: real64
use accelnet_predictor, only: predictor_model,load_predictor_from_networks
use accelnet_descriptors, only: atomic_structure,neighbor_data,build_neighbor_list
use accelnet_batch_target_serial, only: target_model,target_workspace,target_profile,evaluate_batch_target
use batch_test_support, only: make_structure
implicit none
type(predictor_model)::model
type(target_model)::packed
type(target_workspace)::work
type(target_profile)::p
type(atomic_structure)::s
type(neighbor_data)::nb
character(len=1024)::files(2),directory
real(real64)::e(1),sumtimes(7),f(3,192)
integer::a,r,first,last
call get_command_argument(1,directory)
files(1)=trim(directory)//'/Ti.nn.ascii';files(2)=trim(directory)//'/O.nn.ascii'
call load_predictor_from_networks(files,model)
call packed%initialize(model,use_host=.true.)
call make_structure(192,2,s)
call build_neighbor_list(s,model%maximum_cutoff,nb,model%minimum_distance)
sumtimes=0
 do r=1,101
 do a=1,192
 first=nb%offsets(a);last=nb%offsets(a+1)-1;f=0
 call evaluate_batch_target(packed,s%species,[a],[first,last+1],nb%atom_indices,nb%displacements,e,f,work, &
 profile=p,energy_only=.true.)
 if(r>1) sumtimes=sumtimes+[p%prepare,p%upload,p%descriptors,p%network,p%forces,p%download,p%total]
 end do
 end do
print *,sumtimes/100*1000
end program
