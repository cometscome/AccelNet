/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This file serves as an interface to the AccelNet package.

   Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Michael S. Chen, Markland Group, Stanford University
------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include "pair_accelnet.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

extern "C"
{
#include "accelnet.h"
}
  
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairAccelNet::PairAccelNet(LAMMPS *lmp) : Pair(lmp), cut_global(0.0), stat(0),
  chebyshev_mode(ACCELNET_CHEBYSHEV_AUTO), initialized(false),
  atom_types(NULL), pot_files(NULL)
{
  
}

/* ---------------------------------------------------------------------- */

PairAccelNet::~PairAccelNet()
{
  if (initialized) {
    accelnet_final(&stat);
  }

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

  }
  if (atom_types) memory->destroy(atom_types);
  if (pot_files) memory->destroy(pot_files);
}

/* ---------------------------------------------------------------------- */

void PairAccelNet::compute(int eflag, int vflag)
{
  int i,ii,inum,itype;
  int j,jj,jnum;
  int *ilist,*numneigh,**firstneigh;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    double E_i = 0.0;
    i = ilist[ii];
    itype = type[i];
    double icoo[3] = { x[i][0], x[i][1], x[i][2] };
      
    jnum = numneigh[i];
    std::vector<int> jtype(jnum);
    std::vector<int> jlist(jnum);
    std::vector<double> jcoo(3*jnum);
    for (jj = 0; jj < jnum; jj++) {
      j = firstneigh[i][jj];
      j &= NEIGHMASK;
      jlist[jj] = j + 1;
      jtype[jj] = type[j];
      jcoo[3*jj] = x[j][0];
      jcoo[3*jj+1] = x[j][1];
      jcoo[3*jj+2] = x[j][2];
    }

    accelnet_atomic_energy_and_forces(icoo, itype, i+1,
				   jnum, jnum ? &jcoo[0] : NULL,
				   jnum ? &jtype[0] : NULL,
				   jnum ? &jlist[0] : NULL,
				   atom->nmax, &E_i,
				   (double*)&(f[0][0]), &stat);
    
    if (stat != 0) {
      snprintf(error_buffer,sizeof(error_buffer),"AccelNet error code: %d",stat);
      error->all(FLERR,error_buffer);
    }

    if (evflag) ev_tally(0,0,nlocal, 1,
			 E_i,0.0,0.0,0.0,0.0,0.0);
    
  }
  
  if (vflag_fdotr) virial_fdotr_compute();
  
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairAccelNet::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairAccelNet::settings(int narg, char **arg)
{
  int model_offset = 0;
  chebyshev_mode = ACCELNET_CHEBYSHEV_AUTO;
  if (narg == atom->ntypes + 1) {
    model_offset = 1;
    if (strcmp(arg[0],"auto") == 0)
      chebyshev_mode = ACCELNET_CHEBYSHEV_AUTO;
    else if (strcmp(arg[0],"direct") == 0)
      chebyshev_mode = ACCELNET_CHEBYSHEV_DIRECT;
    else if (strcmp(arg[0],"moment") == 0)
      chebyshev_mode = ACCELNET_CHEBYSHEV_MOMENT;
    else
      error->all(FLERR,"AccelNet evaluation mode must be auto, direct, or moment");
  } else if (narg != atom->ntypes) {
    error->all(FLERR,"pair_style accelnet requires [auto|direct|moment] and one potential per atom type");
  }

  memory->create(atom_types, atom->ntypes, 17, "pair:atom_types");
  memory->create(pot_files, atom->ntypes, 1025, "pair:pot_files");

  for (int i = 0; i < atom->ntypes; i++) {
    const char *potential = arg[i + model_offset];
    if (strlen(potential) > 1024)
      error->all(FLERR,"AccelNet potential path is too long");
    snprintf(pot_files[i],1025,"%s",potential);

    std::string filename(potential);
    std::string::size_type slash = filename.find_last_of("/\\");
    std::string basename = slash == std::string::npos ? filename : filename.substr(slash+1);
    std::string::size_type dot = basename.find('.');
    std::string species = basename.substr(0,dot);
    if (species.empty() || species.size() > 16)
      error->all(FLERR,"Cannot infer AccelNet atom type from potential filename");
    snprintf(atom_types[i],17,"%s",species.c_str());
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairAccelNet::coeff(int narg, char **arg)
{
  if (narg != 2)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR, arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR, arg[1],atom->ntypes,jlo,jhi);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAccelNet::init_style()
{
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  
  if (!initialized) {
    accelnet_init(atom->ntypes, atom_types, &stat);
    if (stat == ACCELNET_OK) initialized = true;
  }
  if (stat != 0) {
      snprintf(error_buffer,sizeof(error_buffer),"AccelNet error code: %d",stat);
      error->all(FLERR,error_buffer);
  }
  accelnet_set_chebyshev_evaluation(chebyshev_mode, &stat);
  if (stat != ACCELNET_OK ||
      accelnet_get_chebyshev_evaluation() != chebyshev_mode) {
    snprintf(error_buffer,sizeof(error_buffer),
             "AccelNet failed to select Chebyshev evaluation mode (error code: %d)",stat);
    error->all(FLERR,error_buffer);
  }
  if (!accelnet_all_loaded()) {    
    for (int i = 0; i < atom->ntypes; i++) {
      accelnet_load_potential(i+1, pot_files[i], &stat);
      if (stat != ACCELNET_OK) {
        snprintf(error_buffer,sizeof(error_buffer),
                 "AccelNet failed to load potential %s (error code: %d)",pot_files[i],stat);
        error->all(FLERR,error_buffer);
      }
    }
    if (!accelnet_all_loaded())
      error->all(FLERR,"AccelNet did not load all potentials");
  }

  cut_global = accelnet_Rc_max;
  if (comm->me == 0) {
    const char *mode_name = chebyshev_mode == ACCELNET_CHEBYSHEV_DIRECT ? "direct" :
                            chebyshev_mode == ACCELNET_CHEBYSHEV_MOMENT ? "moment" : "auto";
    if (screen) fprintf(screen,"AccelNet Chebyshev evaluation mode: %s\n",mode_name);
    if (logfile) fprintf(logfile,"AccelNet Chebyshev evaluation mode: %s\n",mode_name);
  }
  
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAccelNet::init_one(int i, int j)
{
  return cut_global;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairAccelNet::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {}
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairAccelNet::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {}
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairAccelNet::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairAccelNet::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */
