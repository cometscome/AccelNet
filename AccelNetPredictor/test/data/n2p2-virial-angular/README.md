# Synthetic angular virial fixture

This hand-authored H/O model tests inference, not physical accuracy. It uses
a linear output layer with nonzero radial and angular weights, nontrivial
descriptor scaling, energy/length normalization, and atomic reference energies.

- H: G2 and G5 (n2p2 SF types 2 and 9), eligible for direct contraction.
- O: G2, G4, and G5 (types 2, 3, and 9), requiring the full-Jacobian path.

Using both species together exercises different derivative paths in the atomic
API and the structure API. The tests assert that those paths are active and
compare energies, forces, and virials. Three-atom molecular and periodic cases
give nonzero angular contributions for the finite-difference tests.

All files are self-contained; no external n2p2 installation is needed.
