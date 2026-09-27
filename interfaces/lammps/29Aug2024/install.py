#!/usr/bin/env python3
"""Install the AccelNet CPU/GPU interface into a LAMMPS 29Aug2024 source tree."""
import argparse
import shutil
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('lammps', type=Path)
a = parser.parse_args()
root = a.lammps.resolve()
here = Path(__file__).resolve().parent
repo = here.parents[2]
if '29 Aug 2024' not in (root/'src/version.h').read_text():
    parser.error('This interface targets LAMMPS 29 Aug 2024 (Update 4); use a matching source tree')
cmake = root / 'cmake/CMakeLists.txt'
text = cmake.read_text()
if '\n  ACCELNET\n' not in text:
    text = text.replace('set(STANDARD_PACKAGES\n', 'set(STANDARD_PACKAGES\n  ACCELNET\n')
if 'foreach(PKG_WITH_INCL ACCELNET ' not in text:
    text = text.replace('foreach(PKG_WITH_INCL KSPACE ', 'foreach(PKG_WITH_INCL ACCELNET KSPACE ')
cmake.write_text(text)
for src, dst in [('ACCELNET', 'src/ACCELNET'), ('GPU', 'src/GPU'), ('lib-gpu', 'lib/gpu')]:
    shutil.copytree(here/src, root/dst, dirs_exist_ok=True, copy_function=shutil.copy)
for name in ['accelnet.h', 'accelnet_target.h']:
    shutil.copy(repo/'AccelNetPredictor/include'/name, root/'src/ACCELNET'/name)
shutil.copy(here/'cmake/ACCELNET.cmake', root/'cmake/Modules/Packages/ACCELNET.cmake')

# GPU auto-binning in this release reads orthogonal sublo/subhi even for
# triclinic domains, where those fields are not initialized. Preserve the
# normal bbox-based sorter for triclinic boxes (including neigh no).
atom_cpp = root/'src/atom.cpp'
text = atom_cpp.read_text()
text = text.replace('#ifdef LMP_GPU\n  if (userbinsize == 0.0) {',
                    '#ifdef LMP_GPU\n  if (userbinsize == 0.0 && !domain->triclinic) {')
atom_cpp.write_text(text)
