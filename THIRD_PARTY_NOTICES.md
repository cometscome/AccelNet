# Third-party notices

Most source code in AccelNet is licensed under the MIT License in `LICENSE`.
The following files contain code derived from other projects and retain their
original licenses and copyright notices.

## ænet linked-cell implementation

File:

- `AccelNetDescriptors/src/accelnet_legacy_lcl.f90`

This file is derived from `src/ext/lclist.f90` in ænet. It remains licensed
under the Mozilla Public License 2.0 and retains the original copyright notice
for Nongnuch Artrith and Alexander Urban. See `LICENSES/MPL-2.0.txt`.

## LAMMPS pair style

Files:

- `interfaces/lammps/4Feb2020/USER-ACCELNET/pair_accelnet.cpp`
- `interfaces/lammps/4Feb2020/USER-ACCELNET/pair_accelnet.h`
- `interfaces/lammps/29Aug2024/ACCELNET/pair_accelnet.cpp`
- `interfaces/lammps/29Aug2024/ACCELNET/pair_accelnet.h`

These files implement a pair style using LAMMPS source interfaces and retain
the LAMMPS copyright and GNU General Public License version 2 notices. See
`LICENSES/GPL-2.0-only.txt`. AccelNet's independent Fortran libraries remain
under their file-specific MIT or MPL-2.0 terms; a distributed LAMMPS executable
built with this pair style is subject to the LAMMPS GPL terms.
