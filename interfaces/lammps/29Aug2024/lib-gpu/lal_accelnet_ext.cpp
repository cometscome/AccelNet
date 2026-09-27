// Adapter for LAMMPS 29Aug2024 lib/gpu. Forces are synchronously returned by
// Fortran in Pair::compute, never queued for FixGPU::post_force.
#include "lal_device.h"
#include <set>
namespace LAMMPS_AL {
extern Device<PRECISION,ACC_PRECISION> global_device;
}
using namespace LAMMPS_AL;
struct AccelNetNeighbors {
  Answer<PRECISION,ACC_PRECISION> answer;
  Neighbor neighbor;
  bool initialized = false;
  ~AccelNetNeighbors() {
    neighbor.clear();
    answer.clear();
    if (initialized) global_device.clear();
  }
};
int accelnet_gpu_device(int &device, int &mode) {
#if !defined(USE_CUDA) || !defined(_DOUBLE_DOUBLE)
  return -4;
#else
  if (!global_device.gpu) return -1;
  if (global_device.particle_split() != 1.0) return -8;
  device = global_device.gpu->device_num();
  // NVHPC borrows the primary context already selected by lib/gpu and caches
  // kernels/pools beyond a LAMMPS clear. Keep one reference per used device for
  // the process lifetime, so FixGPU teardown cannot invalidate that cache.
  // Deliberately no static-destructor release: its ordering relative to the
  // OpenMP runtime's exit handlers is unspecified. CUDA reclaims it on exit.
  static std::set<int> retained_devices;
  if (!retained_devices.count(device)) {
    CUcontext context;
    if (cuDevicePrimaryCtxRetain(&context,device) != CUDA_SUCCESS) return -6;
    retained_devices.insert(device);
  }
  mode = global_device.gpu_mode();
  return 0;
#endif
}
void *accelnet_gpu_neighbors_create(int nlocal, int nall, int maxnb,
                                  double cutoff, int &status) {
  auto *ctx = new AccelNetNeighbors;
  status = global_device.init(ctx->answer,false,false,nlocal,nall,0);
  if (!status) {
    ctx->initialized = true;
    // Force one thread per atom for the documented unpacked matrix layout.
    status = global_device.init_nbor(&ctx->neighbor,nlocal,0,nall,0,0,
                                    maxnb,cutoff,false,1);
  }
  if (status) { delete ctx; return nullptr; }
  return ctx;
}
void accelnet_gpu_neighbors_destroy(void *ptr) {
  delete static_cast<AccelNetNeighbors *>(ptr);
}
int accelnet_gpu_neighbors(void *ptr, int ago, int nlocal, int nall,
                          double **x, int *type, double *lo, double *hi,
                          void *tags, const void *&device_x,
                          const void *&device_nb, int &pitch, int &maxnb) {
  auto &ctx = *static_cast<AccelNetNeighbors *>(ptr);
  auto &atom = global_device.atom;
  bool success = true;
  atom.resize(nall,success);
  ctx.answer.resize(nlocal,success);
  if (!success) return -3;
  if (nlocal == 0) { pitch = 0; maxnb = 0; return 0; }
  atom.data_unavail();
  atom.cast_copy_x(x,type);
  if (ago == 0) {
    ctx.neighbor.resize(nlocal,0,ctx.neighbor.max_nbors(),success);
    if (!success) return -3;
    int mn;
    ctx.neighbor.build_nbor_list(x,nlocal,0,nall,atom,lo,hi,
        static_cast<tagint *>(tags),nullptr,nullptr,success,mn,ctx.answer.error_flag);
    if (!success) return -3;
  }
  global_device.gpu->sync();
  if (ctx.answer.error_flag[0]) return -3;
  device_x = reinterpret_cast<const void *>(atom.x.device.begin());
  device_nb = reinterpret_cast<const void *>(ctx.neighbor.dev_nbor.begin());
  pitch = ctx.neighbor.nbor_pitch();
  maxnb = ctx.neighbor.max_nbors();
  return 0;
}
