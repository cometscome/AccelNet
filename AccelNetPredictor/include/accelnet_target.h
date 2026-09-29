#ifndef ACCELNET_TARGET_H
#define ACCELNET_TARGET_H
/* Explicit CPU execution for tests/host threading. Nonnegative device IDs
   require actual GPU execution and never silently fall back to the CPU. */
#define ACCELNET_TARGET_HOST (-1)
#ifdef __cplusplus
extern "C" {
#endif
/* Independent instance; caller owns the handle and must destroy it before the
   CUDA context is released. Paths may be supplied in any order; species IDs
   always follow embedded model ordering (query them below). FP64 only.
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
/* As above, with Chebyshev convention 0, 1, or 10. The older constructors
   retain version 0. This convention is independent of direct/moment mode. */
int accelnet_target_create_versioned(int nspecies, const char *const *paths, int device,
                                     int chebyshev_version, int chebyshev_mode,
                                     int g5_mode, double *cutoff, void **handle, char *error);
/* Load a supported n2p2 2G model directory directly, without conversion.
   nspecies is returned; query each one-based species ID below for its symbol.
   Invalid arguments, malformed/unsupported model contents, and missing required
   files return nonzero status, a diagnostic, a null handle and zero outputs. */
int accelnet_target_create_n2p2(const char *directory, int device, int g5_mode,
                               int *nspecies, double *cutoff, void **handle, char *error);
/* capacity includes the null terminator. Valid handles from any constructor
   are accepted. A capacity of 17 suffices for the model's element names. */
int accelnet_target_get_species(void *handle, int species, int capacity,
                                char *symbol, char *error);
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
