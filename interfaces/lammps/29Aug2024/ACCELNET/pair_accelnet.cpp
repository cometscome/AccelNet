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
#include "utils.h"

extern "C"
{
#include "accelnet.h"
}
  
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairAccelNet::PairAccelNet(LAMMPS *lmp) : Pair(lmp), cut_global(0.0), stat(0),
  chebyshev_evaluation_mode(ACCELNET_CHEBYSHEV_AUTO),
  g5_evaluation_mode(ACCELNET_G5_AUTO), initialized(false), n2p2_mode(false),
  type_map(NULL), atom_types(NULL), pot_files(NULL), n2p2_directory(NULL)
{
  manybody_flag = 1;
  one_coeff = 1;
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
  if (type_map) memory->destroy(type_map);
  if (n2p2_directory) memory->destroy(n2p2_directory);
}

/* ---------------------------------------------------------------------- */

void PairAccelNet::compute(int eflag, int vflag)
{
  int i,ii,inum,itype;
  int j,jj,jnum;
  int *ilist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  std::vector<int> jtype;
  std::vector<int> jlist;
  std::vector<double> jcoo;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    double E_i = 0.0;
    i = ilist[ii];
    itype = type_map[type[i]];
    double icoo[3] = { x[i][0], x[i][1], x[i][2] };
      
    jnum = numneigh[i];
    jtype.resize(jnum);
    jlist.resize(jnum);
    jcoo.resize(3*jnum);
    for (jj = 0; jj < jnum; jj++) {
      j = firstneigh[i][jj];
      j &= NEIGHMASK;
      jlist[jj] = j + 1;
      jtype[jj] = type_map[type[j]];
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
  const int ntypes = atom->ntypes;
  int cursor = 0;
  chebyshev_evaluation_mode = ACCELNET_CHEBYSHEV_AUTO;
  g5_evaluation_mode = ACCELNET_G5_AUTO;
  n2p2_mode = false;
  if (narg > 0 && (strcmp(arg[0],"auto") == 0 || strcmp(arg[0],"direct") == 0 ||
                   strcmp(arg[0],"moment") == 0)) {
    cursor = 1;
    if (strcmp(arg[0],"auto") == 0)
      chebyshev_evaluation_mode = ACCELNET_CHEBYSHEV_AUTO;
    else if (strcmp(arg[0],"direct") == 0)
      chebyshev_evaluation_mode = ACCELNET_CHEBYSHEV_DIRECT;
    else if (strcmp(arg[0],"moment") == 0)
      chebyshev_evaluation_mode = ACCELNET_CHEBYSHEV_MOMENT;
    else
      error->all(FLERR,"AccelNet Chebyshev mode must be auto, direct, or moment");
  }

  int remaining = narg - cursor;
  n2p2_mode = remaining >= 1 && strcmp(arg[cursor],"n2p2") == 0;
  const int base_arguments = n2p2_mode ? ntypes + 2 : ntypes;
  if (remaining != base_arguments && remaining != base_arguments + 2)
    error->all(FLERR,"Expected [mode] potentials... [g5 MODE] or [mode] n2p2 directory elements... [g5 MODE]");
  if (remaining == base_arguments + 2) {
    if (strcmp(arg[cursor + base_arguments],"g5") != 0)
      error->all(FLERR,"Expected 'g5 MODE' after the AccelNet model arguments");
    const char *g5_mode = arg[cursor + base_arguments + 1];
    if (strcmp(g5_mode,"auto") == 0)
      g5_evaluation_mode = ACCELNET_G5_AUTO;
    else if (strcmp(g5_mode,"direct") == 0)
      g5_evaluation_mode = ACCELNET_G5_DIRECT;
    else if (strcmp(g5_mode,"moment") == 0)
      g5_evaluation_mode = ACCELNET_G5_MOMENT_FORCE;
    else
      error->all(FLERR,"AccelNet g5 mode must be auto, direct, or moment");
  }

  memory->create(atom_types, atom->ntypes, 17, "pair:atom_types");
  memory->create(type_map, atom->ntypes + 1, "pair:type_map");
  type_map[0] = 0;

  if (n2p2_mode) {
    const char *directory = arg[cursor + 1];
    if (strlen(directory) > 1024)
      error->all(FLERR,"AccelNet n2p2 model directory path is too long");
    memory->create(n2p2_directory, 1025, "pair:n2p2_directory");
    snprintf(n2p2_directory,1025,"%s",directory);
    for (int i = 0; i < ntypes; i++) {
      const char *species = arg[cursor + 2 + i];
      if (strlen(species) == 0 || strlen(species) > 16)
        error->all(FLERR,"Invalid element name in AccelNet n2p2 mapping");
      snprintf(atom_types[i],17,"%s",species);
    }
    return;
  }

  memory->create(pot_files, atom->ntypes, 1025, "pair:pot_files");

  for (int i = 0; i < ntypes; i++) {
    const char *potential = arg[cursor + i];
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
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);
  
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
  neighbor->add_request(this, NeighConst::REQ_FULL);
  
  if (!initialized) {
    if (n2p2_mode)
      accelnet_init_n2p2(n2p2_directory, &stat);
    else
      accelnet_init(atom->ntypes, atom_types, &stat);
    if (stat == ACCELNET_OK) initialized = true;
  }
  if (stat != 0) {
      snprintf(error_buffer,sizeof(error_buffer),"AccelNet error code: %d",stat);
      error->all(FLERR,error_buffer);
  }
  if (n2p2_mode) {
    std::vector<int> input_types(atom->ntypes);
    for (int i = 0; i < atom->ntypes; i++) input_types[i] = i + 1;
    accelnet_convert_atom_types(atom->ntypes, atom_types, atom->ntypes,
                                input_types.data(), type_map + 1, &stat);
    if (stat != ACCELNET_OK)
      error->all(FLERR,"LAMMPS atom types do not match elements in the n2p2 model");
  } else {
    for (int i = 1; i <= atom->ntypes; i++) type_map[i] = i;
  }
  if (!n2p2_mode && !accelnet_all_loaded()) {
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

  accelnet_set_chebyshev_evaluation(chebyshev_evaluation_mode, &stat);
  if (stat != ACCELNET_OK) {
    snprintf(error_buffer,sizeof(error_buffer),
             "AccelNet failed to set Chebyshev evaluation mode (error code: %d)",stat);
    error->all(FLERR,error_buffer);
  }
  accelnet_set_g5_evaluation(g5_evaluation_mode, &stat);
  if (stat != ACCELNET_OK) {
    snprintf(error_buffer,sizeof(error_buffer),
             "AccelNet failed to set G5 evaluation mode (error code: %d)",stat);
    error->all(FLERR,error_buffer);
  }
  if (comm->me == 0) {
    const char *chebyshev_mode_name = chebyshev_evaluation_mode == ACCELNET_CHEBYSHEV_DIRECT ? "direct" :
                                      chebyshev_evaluation_mode == ACCELNET_CHEBYSHEV_MOMENT ? "moment" : "auto";
    const char *mode_name = g5_evaluation_mode == ACCELNET_G5_DIRECT ? "direct" :
                            g5_evaluation_mode == ACCELNET_G5_MOMENT_FORCE ? "moment" : "auto";
    if (screen) fprintf(screen,"AccelNet Chebyshev evaluation mode: %s\n",chebyshev_mode_name);
    if (logfile) fprintf(logfile,"AccelNet Chebyshev evaluation mode: %s\n",chebyshev_mode_name);
    if (screen) fprintf(screen,"AccelNet G5 evaluation mode: %s\n",mode_name);
    if (logfile) fprintf(logfile,"AccelNet G5 evaluation mode: %s\n",mode_name);
    if (n2p2_mode) {
      if (screen) fprintf(screen,"AccelNet n2p2 model directory: %s\n",n2p2_directory);
      if (logfile) fprintf(logfile,"AccelNet n2p2 model directory: %s\n",n2p2_directory);
    }
  }

  cut_global = accelnet_Rc_max;
  
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
