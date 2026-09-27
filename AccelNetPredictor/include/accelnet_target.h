#ifndef ACCELNET_TARGET_H
#define ACCELNET_TARGET_H
#ifdef __cplusplus
extern "C" {
#endif
/* Independent instance; caller owns the handle and must destroy it before the
   CUDA context is released. Paths follow embedded species ordering. FP64 only.
   CSR indices/offsets are one-based. Forces and column-major virial are additive.
   Error buffer must have 512 bytes. Status zero means success. */
int accelnet_target_create(int nspecies, const char *const *paths, int device,
                           int mode, double *cutoff, void **handle, char *error);
/* Independent modes: Chebyshev 0 auto / 1 direct / 2 moment;
   G5 0 auto (16 angular neighbors) / 1 direct / 2 moment with 16-neighbor threshold /
   3 forced moment. Ineligible fractional/high G5 powers remain direct. */
int accelnet_target_create_modes(int nspecies, const char *const *paths, int device,
                                 int chebyshev_mode, int g5_mode, double *cutoff,
                                 void **handle, char *error);
void accelnet_target_destroy(void *handle);
int accelnet_target_compute(void *handle, int natoms, int nrows, int nedges,
                            const int *species, const int *centers,
                            const int *offsets, const int *indices,
                            const double *dr, double *energies, double *forces,
                            double *virial, char *error);
/* Borrowed lib/gpu CUDA arrays: double4 x and unpacked neighbor matrix,
   threads_per_atom=1, pitch entries per row. Caller synchronizes producer.
   Supported only for validated atomic models/types and full neighbor lists. */
int accelnet_target_compute_device(void *handle, int natoms, int nrows,
                                   int maxnb, int pitch, const void *x,
                                   const void *neighbors, double *energies,
                                   double *forces, double *virial);
#ifdef __cplusplus
}
#endif
#endif
