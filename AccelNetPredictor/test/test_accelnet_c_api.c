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
        accelnet_final(&stat);
        if (!check(stat == ACCELNET_OK, "C n2p2 final")) return 1;
        accelnet_init_n2p2(argv[3], &stat);
        if (!check(stat == ACCELNET_OK && accelnet_all_loaded(), "C one-shot n2p2 init")) return 1;
        accelnet_final(&stat);
        if (!check(stat == ACCELNET_OK, "C one-shot n2p2 final")) return 1;
    }
    return 0;
}
