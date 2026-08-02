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
  ~PairAccelNet() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

 protected:
  void allocate();
  double cut_global;

  int stat;
  int g5_evaluation_mode;
  bool initialized;
  char **atom_types;
  char **pot_files;
  char error_buffer[128];
  
};

}

#endif
#endif
