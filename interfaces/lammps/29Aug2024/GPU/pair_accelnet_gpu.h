#ifdef PAIR_CLASS
PairStyle(accelnet/gpu,PairAccelNetGPU)
#else
#ifndef LMP_PAIR_ACCELNET_GPU_H
#define LMP_PAIR_ACCELNET_GPU_H
#include "pair_accelnet.h"
#include <vector>
namespace LAMMPS_NS {
class PairAccelNetGPU : public PairAccelNet {
 public:
  PairAccelNetGPU(class LAMMPS *);
  ~PairAccelNetGPU() override;
  void settings(int,char **) override;
  void init_style() override;
  void compute(int,int) override;
 private:
  void *target_context = nullptr, *neighbor_context = nullptr;
  int gpu_mode = 0;
  std::vector<int> species,centers,offsets,indices;
  std::vector<double> displacements,energies;
};
}
#endif
#endif
