/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef PAIR_CLASS

PairStyle(accelnet,PairAccelNet)

#else

#ifndef LMP_PAIR_ACCELNET
#define LMP_PAIR_ACCELNET

#include "pair.h"

namespace LAMMPS_NS {

class PairAccelNet : public Pair {
 public:
  PairAccelNet(class LAMMPS *);
  virtual ~PairAccelNet();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 protected:
  virtual void allocate();
  double cut_global;

  int stat;
  int chebyshev_mode;
  bool initialized;
  bool n2p2_mode;
  int *type_map;
  char **atom_types;
  char **pot_files;
  char *n2p2_directory;
  char error_buffer[128];
  
};

}

#endif
#endif
