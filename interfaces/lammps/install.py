#!/usr/bin/env python3
"""Install AccelNet CPU/GPU support into LAMMPS 22Jul2025 Update 6 (recommended)
or 29Aug2024 Update 4. The same adapter sources serve both releases.
"""
import argparse
import re
import shutil
from pathlib import Path

SUPPORTED = {("22 Jul 2025", "Update 6"), ("29 Aug 2024", "Update 4")}


def register_package(text):
    """Check and update both CMake registration sites before writing files."""
    patterns = (
        r"set\(STANDARD_PACKAGES\b[^)]*\)",
        r"foreach\(PKG_WITH_INCL\s+[^)]*\bKSPACE\b[^)]*\)",
    )
    for pattern in patterns:
        matches = list(re.finditer(pattern, text))
        if len(matches) != 1:
            raise ValueError("Unrecognized LAMMPS CMake package registration; no files changed")
        match = matches[0]
        block = match.group()
        if "ACCELNET" not in block.split():
            if block.startswith("set("):
                block = block.replace("STANDARD_PACKAGES", "STANDARD_PACKAGES\n  ACCELNET", 1)
            else:
                block = block.replace("PKG_WITH_INCL", "PKG_WITH_INCL ACCELNET", 1)
            text = text[:match.start()] + block + text[match.end():]
    return text


def install(root):
    here = Path(__file__).resolve().parent
    adapter = here / "29Aug2024"
    repo = here.parents[1]
    version = (root / "src/version.h").read_text()
    release = re.search(r'^#define\s+LAMMPS_VERSION\s+"([^"]+)"', version, re.M)
    update = re.search(r'^#define\s+LAMMPS_UPDATE\s+"([^"]+)"', version, re.M)
    key = (release.group(1) if release else "", update.group(1) if update else "")
    if key not in SUPPORTED:
        raise ValueError(
            "Supported LAMMPS versions: 22 Jul 2025 Update 6 (recommended; "
            "tag stable_22Jul2025_update6) and 29 Aug 2024 Update 4. "
            "Other versions, including 2 Sep 2026, are not supported by this installer. "
            "No files changed."
        )
    cmake = root / "cmake/CMakeLists.txt"
    cmake_text = register_package(cmake.read_text())
    atom_cpp = root / "src/atom.cpp"
    atom_text = atom_cpp.read_text()
    old = "#ifdef LMP_GPU\n  if (userbinsize == 0.0) {"
    new = "#ifdef LMP_GPU\n  if (userbinsize == 0.0 && !domain->triclinic) {"
    if old not in atom_text and new not in atom_text:
        raise ValueError("Unrecognized LAMMPS GPU atom-sorting code; no files changed")
    atom_text = atom_text.replace(old, new)

    # Both releases use the same pair styles and Fortran numerical kernels.
    for src, dst in (("ACCELNET", "src/ACCELNET"), ("GPU", "src/GPU"),
                     ("lib-gpu", "lib/gpu")):
        shutil.copytree(adapter / src, root / dst, dirs_exist_ok=True,
                        copy_function=shutil.copy)
    for name in ("accelnet.h", "accelnet_target.h"):
        shutil.copy(repo / "AccelNetPredictor/include" / name,
                    root / "src/ACCELNET" / name)
    shutil.copy(adapter / "cmake/ACCELNET.cmake",
                root / "cmake/Modules/Packages/ACCELNET.cmake")
    cmake.write_text(cmake_text)
    # Keep the normal bbox sorter for triclinic boxes, including neigh no.
    atom_cpp.write_text(atom_text)
    return " ".join(key)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lammps", type=Path, help="LAMMPS source directory")
    args = parser.parse_args()
    try:
        version = install(args.lammps.resolve())
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Installed AccelNet CPU/GPU interface for LAMMPS {version}")
    print("Build CPU with PKG_ACCELNET=ON, PKG_GPU=OFF; "
          "for GPU use the NVHPC/CUDA recipe in interfaces/lammps/README.md.")


if __name__ == "__main__":
    main()
