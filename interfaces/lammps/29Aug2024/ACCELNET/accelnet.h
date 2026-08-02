#ifndef ACCELNET_H_INCLUDED
#define ACCELNET_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
#  define ACCELNET_BOOL bool
#else
#  include <stdbool.h>
#  define ACCELNET_BOOL _Bool
#endif

void accelnet_init(int ntypes, char *atom_types[], int *stat);
void accelnet_final(int *stat);
void accelnet_print_info(void);
void accelnet_set_chebyshev_version(int version, int *stat);
void accelnet_set_g5_evaluation(int mode, int *stat);
void accelnet_load_potential(int type_id, char *filename, int *stat);
void accelnet_load_potential_ascii(int type_id, char *filename, int *stat);
ACCELNET_BOOL accelnet_all_loaded(void);

double accelnet_free_atom_energy(int type_id);

void accelnet_atomic_energy(const double coo_i[3], int type_i, int n_j,
                            const double coo_j[], const int type_j[],
                            double *energy_i, int *stat);
void accelnet_atomic_energy_and_forces(
    const double coo_i[3], int type_i, int index_i, int n_j,
    const double coo_j[], const int type_j[], const int index_j[],
    int natoms, double *energy_i, double forces[], int *stat);

void accelnet_convert_atom_types(int ntypes_in, char *atom_types[],
                                 int natoms_in, const int type_id_in[],
                                 int type_id_out[], int *stat);

void accelnet_nbl_init(const double lattice[9], int natoms,
                       const int atom_types[], double coordinates[],
                       ACCELNET_BOOL cartesian, ACCELNET_BOOL pbc);
void accelnet_nbl_final(void);
void accelnet_nbl_neighbors(int iatom, int *nnb, double nbcoo[],
                            double nbdist[], int nblist[], int nbtype[]);

void accelnet_sfb_init(int ntypes, char *atom_types[], int radial_order,
                       int angular_order, double radial_cutoff,
                       double angular_cutoff, int *stat);
void accelnet_sfb_final(int *stat);
int accelnet_sfb_nvalues(void);
void accelnet_sfb_eval(int type_i, const double coo_i[3], int n_j,
                       const int type_j[], const double coo_j[], int nvalues,
                       double values[], int *stat);
void accelnet_sfb_reconstruct_radial(int nvalues, const double values[],
                                     int nx, double x[], double y[], int *stat);

extern int ACCELNET_OK;
extern int ACCELNET_ERR_INIT;
extern int ACCELNET_ERR_MALLOC;
extern int ACCELNET_ERR_IO;
extern int ACCELNET_ERR_TYPE;
extern int ACCELNET_ERR_ARGUMENT;
extern int ACCELNET_TYPELEN;
extern int ACCELNET_PATHLEN;
extern ACCELNET_BOOL ACCELNET_TRUE;
extern ACCELNET_BOOL ACCELNET_FALSE;

extern int ACCELNET_G5_AUTO;
extern int ACCELNET_G5_DIRECT;
extern int ACCELNET_G5_MOMENT;
extern int ACCELNET_G5_MOMENT_FORCE;

extern int accelnet_nsf_max;
extern int accelnet_nnb_max;
extern double accelnet_Rc_min;
extern double accelnet_Rc_max;

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* ACCELNET_H_INCLUDED */
