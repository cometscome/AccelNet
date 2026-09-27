#include "accelnet.h"

#include <math.h>
#include <stdio.h>

static int check(int condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "FAILED: %s\n", message);
        return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    char *species[] = {"Ti", "O"};
    double center[3] = {0.0, 0.0, 0.0};
    double neighbors[3] = {1.9, 0.0, 0.0};
    int neighbor_types[1] = {2};
    int stat = -1;
    double energy = 0.0;

    if (argc != 4) return 2;
    accelnet_init(2, species, &stat);
    if (!check(stat == ACCELNET_OK, "C init")) return 1;
    if (!check(accelnet_get_chebyshev_evaluation() == ACCELNET_CHEBYSHEV_AUTO,
               "C Chebyshev default")) return 1;
    accelnet_set_chebyshev_evaluation(ACCELNET_CHEBYSHEV_DIRECT, &stat);
    if (!check(stat == ACCELNET_OK &&
               accelnet_get_chebyshev_evaluation() == ACCELNET_CHEBYSHEV_DIRECT,
               "C Chebyshev direct")) return 1;
    accelnet_set_chebyshev_evaluation(3, &stat);
    if (!check(stat == ACCELNET_ERR_ARGUMENT, "C invalid Chebyshev mode")) return 1;
    accelnet_set_chebyshev_evaluation(ACCELNET_CHEBYSHEV_MOMENT, &stat);
    if (!check(stat == ACCELNET_OK, "C Chebyshev moment")) return 1;
    accelnet_load_potential_ascii(1, argv[1], &stat);
    if (!check(stat == ACCELNET_OK, "C load Ti")) return 1;
    accelnet_load_potential_ascii(2, argv[2], &stat);
    if (!check(stat == ACCELNET_OK && accelnet_all_loaded(), "C load O")) return 1;
    if (!check(accelnet_get_chebyshev_evaluation() == ACCELNET_CHEBYSHEV_MOMENT,
               "C Chebyshev mode retained after load")) return 1;
    accelnet_atomic_energy(center, 1, 1, neighbors, neighbor_types, &energy, &stat);
    if (!check(stat == ACCELNET_OK && isfinite(energy), "C atomic energy")) return 1;
    {
        int types[2] = {1, 2}, centers[2] = {1, 2}, offsets[3] = {1, 2, 3};
        int indices[2] = {2, 1};
        double dr[6] = {1.9, 0, 0, -1.9, 0, 0}, energies[2], reference[2];
        double f[6] = {1,2,3,4,5,6}, expected[6] = {1,2,3,4,5,6};
        for (int row = 0; row < 2; ++row) {
            int target_type = types[1-row];
            accelnet_atomic_energy_and_forces(center, types[row], centers[row], 1,
                dr+3*row, &target_type, indices+row, 2, reference+row, expected, &stat);
            if (!check(stat == ACCELNET_OK, "atomic force reference")) return 1;
        }
        accelnet_batch_energy_and_forces(2,2,2,types,centers,offsets,indices,dr,energies,f,&stat);
        if (!check(stat == ACCELNET_OK, "CSR batch status")) return 1;
        for (int i=0;i<2;++i)
            if (!check(fabs(energies[i]-reference[i])<1e-10, "CSR batch energies")) return 1;
        for (int i=0;i<6;++i)
            if (!check(fabs(f[i]-expected[i])<1e-10, "CSR additive forces")) return 1;
        offsets[1]=4;
        accelnet_batch_energy_and_forces(2,2,2,types,centers,offsets,indices,dr,energies,f,&stat);
        if (!check(stat == ACCELNET_ERR_ARGUMENT, "reject invalid CSR")) return 1;
        int empty_offset=1;
        accelnet_batch_energy_and_forces(0,0,0,NULL,NULL,&empty_offset,NULL,NULL,NULL,NULL,&stat);
        if (!check(stat == ACCELNET_OK, "empty CPU batch")) return 1;
    }
    if (!check(accelnet_nsf_max > 0 && accelnet_Rc_max > 0.0, "C globals")) return 1;
    accelnet_final(&stat);
    if (!check(stat == ACCELNET_OK && !accelnet_all_loaded(), "C final")) return 1;

    {
        char *hydrogen[] = {"H"};
        double n2p2_center[3] = {0.0, 0.0, 0.0};
        double n2p2_neighbor[3] = {2.0, 0.0, 0.0};
        int n2p2_type[1] = {1};
        accelnet_init(1, hydrogen, &stat);
        if (!check(stat == ACCELNET_OK, "C n2p2 init")) return 1;
        accelnet_load_n2p2(argv[3], &stat);
        if (!check(stat == ACCELNET_OK && accelnet_all_loaded(), "C n2p2 load")) return 1;
        accelnet_atomic_energy(n2p2_center, 1, 1, n2p2_neighbor, n2p2_type, &energy, &stat);
        if (!check(stat == ACCELNET_OK && isfinite(energy), "C n2p2 atomic energy")) return 1;
        {
            int types[2] = {1, 1}, center_id = 1, offsets[2] = {1, 2}, index = 2;
            double batch_energy, forces[6] = {0};
            accelnet_batch_energy_and_forces(2,1,1,types,&center_id,offsets,&index,
                n2p2_neighbor,&batch_energy,forces,&stat);
            if (!check(stat == ACCELNET_OK && fabs(batch_energy-energy)<1e-10,
                       "CSR metadata refreshed after model replacement")) return 1;
        }
        accelnet_final(&stat);
        if (!check(stat == ACCELNET_OK, "C n2p2 final")) return 1;
        accelnet_init_n2p2(argv[3], &stat);
        if (!check(stat == ACCELNET_OK && accelnet_all_loaded(), "C one-shot n2p2 init")) return 1;
        accelnet_final(&stat);
        if (!check(stat == ACCELNET_OK, "C one-shot n2p2 final")) return 1;
    }
    return 0;
}
