# Changes

## 1.1.0 — 2026-09-28

### Shared CPU/GPU inference

- Production structure/file, atomic Fortran/C, CSR batch, ænet-compatible SFB
  and LAMMPS inference now use the same maintained Fortran numerical kernels.
  Independent legacy evaluators remain available for reference tests and
  standalone low-level descriptor/NN compatibility APIs.
- Serial CPU remains the default, with OpenMP compiled out. The optional target
  API supports explicit CPU OpenMP threading and NVIDIA GPU offload. Ordinary
  CLI/API calls and LAMMPS `accelnet` remain serial within each caller/MPI rank.
- GPU evaluation includes descriptors, neural networks, force contraction and
  virials in FP64. Persistent workspaces reduce repeated allocation and model
  transfers. Host force accumulation avoids fine-grained atomic updates.

### Descriptors and models

- Support ænet Chebyshev, Behler G1--G5 and the AccelNet LJ extension through
  the common backend, including multiple Chebyshev blocks per element and
  mixed descriptor families in multi-element models.
- Extend n2p2 2G-HDNNP support to types **2/3/9/12/13/20--25**, including
  per-element network topology, activation/scaling conventions and
  `normalize_nodes`.
- Provide exact direct/moment Chebyshev evaluation and exact G5 moments for
  integer powers up to 16 when explicitly selected. G5 auto retains its
  order-10 eligibility limit and neighbor threshold. Other angular functions
  use exact direct evaluation where no finite moment factorization is available.
- See [implementation status](docs/implementation-status.md) and
  [model compatibility](docs/model-compatibility.md) for the supported scope
  and remaining limitations.

### Performance

- Optimize common kernels using differentiated Clenshaw and Horner recurrences,
  shared radial/angular work, grouped descriptor evaluation, fused G4
  value/Jacobian traversal and improved moment construction/contraction.
- Remove redundant atomic-energy environment preparation while retaining the
  common numerical pipeline. In the recorded single-core Ti/O comparison,
  the 512-atom slowdown versus the original CPU implementation fell from about
  11% to 1.4%; 64/192-atom cases were about 2% faster in that run. Both CPU
  binaries compiled OpenMP out and used identical drivers.
- Add performance regression checks and retain numerical comparisons against
  independent implementations and finite differences. Timings are specific
  to their recorded models, hardware and compilers; no universal speedup is
  implied. Equations and optimization history are in
  [speedupmethods.md](speedupmethods.md), document revision 1.14.

### LAMMPS and build documentation

- Recommend **LAMMPS 22 Jul 2025 Update 6**, pinned to
  `stable_22Jul2025_update6`, for CPU and GPU use.
- Add `interfaces/lammps/install.py` to install the common adapter, canonical
  headers, CMake registration and triclinic sorting fix automatically. Retain
  29 Aug 2024 Update 4 compatibility and the old installer entry point.
- Support `pair_style accelnet/gpu` through the CUDA GPU package and Fortran
  OpenMP target. The validated compiler stack is NVHPC 25.3 with CUDA 12.8;
  LAMMPS GPU builds require both `nvfortran` and `nvc++`.
- Keep the traditional CPU workflow first in the README and separate CPU
  OpenMP, GPU toolchain requirements and LAMMPS CPU/GPU instructions. See the
  [LAMMPS guide](interfaces/lammps/README.md) for complete commands.
- The 2 Sep 2026 LAMMPS release candidate requires GPU API changes and is not
  accepted by the installer. AMD/Intel GPU execution remains unvalidated.

### Validation records

- [Common energy kernels, methods 1.13](docs/validation/energy-common-2026-09-27/README.md):
  52 CPU tests, 22 bounds/runtime checks, 35 GPU tests, 21 LAMMPS comparisons,
  H100 memory checks and host 2/8-thread checks at revision `f655fb0`.
- [Atomic-energy preparation, methods 1.14](docs/validation/atomic-energy-preparation-2026-09-28/README.md):
  54 CPU tests, 22 bounds/runtime checks, 35 H100 tests and four performance
  gates passed. That campaign did not rerun LAMMPS.
- [Stable LAMMPS integration and installer](docs/validation/lammps-current-2026-09-28/README.md):
  five installer tests, nine CPU MPI comparisons and 21 H100 integration
  comparisons passed after installation and rebuilding. Cases include
  orthogonal/triclinic cells, empty ranks and CPU/GPU neighbor-list modes.

Each report identifies its source revisions, fixtures, tolerances and build
configuration. Test counts depend on the available external fixtures and do
not describe every fresh checkout. Historical Blackwell results use the
retained 2024 LAMMPS release; the 2025 release was checked on H100.
