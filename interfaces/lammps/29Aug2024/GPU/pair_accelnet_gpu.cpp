#include "pair_accelnet_gpu.h"
#include "accelnet_target.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"
#include <climits>
#include <cstring>

using namespace LAMMPS_NS;
int accelnet_gpu_device(int &,int &);
void *accelnet_gpu_neighbors_create(int,int,int,double,int &);
void accelnet_gpu_neighbors_destroy(void *);
int accelnet_gpu_neighbors(void *,int,int,int,double **,int *,double *,double *,void *,
                          const void *&,const void *&,int &,int &);
PairAccelNetGPU::PairAccelNetGPU(LAMMPS *lmp) : PairAccelNet(lmp) {
  suffix_flag |= Suffix::GPU;
  respa_enable = 0;
  restartinfo = 0;
  single_enable = 0;
  no_virial_fdotr_compute = 1;
  GPU_EXTRA::gpu_ready(lmp->modify,lmp->error);
}
PairAccelNetGPU::~PairAccelNetGPU() {
  accelnet_target_destroy(target_context);
  accelnet_gpu_neighbors_destroy(neighbor_context);
}
void PairAccelNetGPU::settings(int narg,char **arg) {
  if (target_context) error->all(FLERR,"Recreate pair_style accelnet/gpu to load another model");
  PairAccelNet::settings(narg,arg);
  if (n2p2_mode) error->all(FLERR,"accelnet/gpu requires embedded networks; convert n2p2 with accelnet-model-converter-fortran first");
}
void PairAccelNetGPU::init_style() {
  if (!force->newton_pair) error->all(FLERR,"accelnet/gpu requires newton on");
  if (atom->molecular != Atom::ATOMIC)
    error->all(FLERR,"accelnet/gpu currently requires atom_style atomic");
  if (std::strcmp(update->integrate_style,"verlet") != 0)
    error->all(FLERR,"accelnet/gpu currently requires run_style verlet");
  if (force->pair != this) error->all(FLERR,"accelnet/gpu does not support hybrid pair styles");
  if (neighbor->nex_type || neighbor->nex_group || neighbor->nex_mol || neighbor->includegroup)
    error->all(FLERR,"accelnet/gpu does not support neighbor exclusions or include groups");
  int device;
  int status = accelnet_gpu_device(device,gpu_mode);
  GPU_EXTRA::check_flag(status,error,world);
  if (!target_context) {
    char detail[512] = {};
    status = accelnet_target_create_modes(atom->ntypes,const_cast<const char **>(pot_files),
        device,chebyshev_evaluation_mode,g5_evaluation_mode,&cut_global,&target_context,detail);
    if (status) error->one(FLERR,std::string("accelnet/gpu: ")+detail);
  }
  if (gpu_mode == 0) {
    neighbor->add_request(this,NeighConst::REQ_FULL);
  }
  {
    // lib/gpu also supplies atom-sorting bin sizes in neigh no mode.
    accelnet_gpu_neighbors_destroy(neighbor_context);
    neighbor_context = accelnet_gpu_neighbors_create(atom->nlocal,atom->nlocal+atom->nghost,
        std::max(1,neighbor->oneatom/20),cut_global+neighbor->skin,status);
    GPU_EXTRA::check_flag(status,error,world);
  }
  if (comm->me == 0 && screen)
    fprintf(screen,"AccelNet GPU: CUDA device %d, FP64 OpenMP target, %s neighbors\n",
            device,gpu_mode == 0 ? "CPU" : "GPU package");
}
void PairAccelNetGPU::compute(int eflag,int vflag) {
  ev_init(eflag,vflag);
  if (vflag_atom || cvflag_atom)
    error->all(FLERR,"accelnet/gpu does not yet support per-atom stress");
  const int nall = atom->nlocal+atom->nghost;
  const int nrows = gpu_mode == 0 ? list->inum : atom->nlocal;
  energies.resize(nrows);
  double w[9] = {};
  int status = 0;
  char detail[512] = {};
  if (gpu_mode == 0) {
    species.assign(atom->type,atom->type+nall);
    centers.resize(nrows);
    offsets.resize(nrows+1);
    indices.clear(); displacements.clear();
    offsets[0] = 1;
    for (int row=0;row<nrows;++row) {
      int i = list->ilist[row]; centers[row] = i+1;
      if (indices.size() > static_cast<size_t>(INT_MAX-list->numneigh[i]-1))
        error->one(FLERR,"accelnet/gpu neighbor list exceeds 32-bit indexing");
      for (int k=0;k<list->numneigh[i];++k) {
        int j = list->firstneigh[i][k] & NEIGHMASK;
        indices.push_back(j+1);
        for (int d=0;d<3;++d) displacements.push_back(atom->x[j][d]-atom->x[i][d]);
      }
      offsets[row+1] = indices.size()+1;
    }
    status = accelnet_target_compute(target_context,nall,nrows,indices.size(),species.data(),
        centers.data(),offsets.data(),indices.data(),displacements.data(),energies.data(),
        atom->f[0],w,detail);
  } else {
    double lo[3],hi[3];
    if (domain->triclinic) domain->bbox(domain->sublo_lamda,domain->subhi_lamda,lo,hi);
    else for (int d=0;d<3;++d) { lo[d]=domain->sublo[d]; hi[d]=domain->subhi[d]; }
    const void *x=nullptr,*nb=nullptr;
    int pitch=0,maxnb=0;
    status = accelnet_gpu_neighbors(neighbor_context,neighbor->ago,nrows,nall,atom->x,
        atom->type,lo,hi,atom->tag,x,nb,pitch,maxnb);
    if (status) error->one(FLERR,"accelnet/gpu GPU neighbor construction failed");
    status = accelnet_target_compute_device(target_context,nall,nrows,maxnb,pitch,x,nb,
        energies.data(),atom->f[0],w);
  }
  if (status) error->one(FLERR,std::string("accelnet/gpu computation failed: ")+detail);
  for (int row=0;row<nrows;++row) {
    int i = gpu_mode == 0 ? centers[row]-1 : row;
    if (eflag_global) eng_vdwl += energies[row];
    if (eflag_atom) eatom[i] += energies[row];
  }
  if (vflag_global) {
    virial[0]+=w[0]; virial[1]+=w[4]; virial[2]+=w[8];
    virial[3]+=w[3]; virial[4]+=w[6]; virial[5]+=w[7];
  }
  // All local and ghost forces are ready here, before Verlet reverse_comm.
}
