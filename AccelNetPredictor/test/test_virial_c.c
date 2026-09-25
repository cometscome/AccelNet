#include "accelnet.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static int check(int condition, const char *label) {
    if (!condition) fprintf(stderr, "FAILED: %s\n", label);
    return condition;
}

static int check_invalid_calls(int initialized) {
    const double center[3] = {0.0, 0.0, 0.0}, neighbors[3] = {0.3, 0.1, 0.2};
    /* center type, center ID, neighbor count, atom count, neighbor type, neighbor ID */
    const int cases[][6] = {
        {1, 0, 1, 2, 1, 2}, {1, 3, 1, 2, 1, 2},
        {1, 1, 1, 2, 1, 0}, {1, 1, 1, 2, 1, 3},
        {0, 1, 1, 2, 1, 2}, {2, 1, 1, 2, 1, 2},
        {1, 1, 1, 2, 0, 2}, {1, 1, 1, 2, 2, 2},
        {1, 1, -1, 2, 1, 2}, {1, 1, 1, 0, 1, 2}
    };
    for (unsigned i=0; i<sizeof cases/sizeof cases[0]; ++i) {
        double forces[6], saved_forces[6], virial[9], saved_virial[9], energy = 123.0;
        int stat, expected = (i >= 4 && i <= 7) ? ACCELNET_ERR_TYPE : ACCELNET_ERR_ARGUMENT;
        if (!initialized) expected = ACCELNET_ERR_INIT;
        for (int j=0; j<6; ++j) forces[j]=saved_forces[j]=0.125+j;
        for (int j=0; j<9; ++j) virial[j]=saved_virial[j]=-0.25+j;
        accelnet_atomic_energy_and_forces_virial(center, cases[i][0], cases[i][1], cases[i][2],
            neighbors, &cases[i][4], &cases[i][5], cases[i][3], &energy, forces, virial, &stat);
        if (!check(stat == expected && energy == 0.0, "C invalid argument/init status and energy")) return 0;
        if (!check(memcmp(forces,saved_forces,sizeof forces) == 0 &&
                   memcmp(virial,saved_virial,sizeof virial) == 0, "C errors preserve both accumulators")) return 0;
    }
    return 1;
}

static int check_cutoff(void) {
    const double center[3] = {0.0, 0.0, 0.0};
    const double factors[] = {0.95, 1.0, 1.05};
    const int type = 1, index = 2;
    double neighbor[3] = {0.0, 0.0, 0.0}, isolated_energy;
    int stat;
    accelnet_atomic_energy(center, 1, 0, neighbor, &type, &isolated_energy, &stat);
    if (!check(stat == ACCELNET_OK, "C isolated reference")) return 0;
    for (int i=0; i<3; ++i) {
        double forces[6] = {0}, virial[9] = {0}, energy, ep, em;
        double r=accelnet_Rc_max*factors[i], h=1e-5;
        neighbor[0]=r;
        accelnet_atomic_energy_and_forces_virial(center, 1, 1, 1, neighbor, &type, &index,
                                                2, &energy, forces, virial, &stat);
        if (!check(stat == ACCELNET_OK && isfinite(energy), "C cutoff evaluation")) return 0;
        neighbor[0]=r*(1+h);
        accelnet_atomic_energy(center, 1, 1, neighbor, &type, &ep, &stat);
        if (!check(stat == ACCELNET_OK, "C cutoff plus")) return 0;
        neighbor[0]=r*(1-h);
        accelnet_atomic_energy(center, 1, 1, neighbor, &type, &em, &stat);
        if (!check(stat == ACCELNET_OK && fabs(virial[0]+(ep-em)/(2*h)) < 1e-7,
                   "C cutoff strain derivative")) return 0;
        if (i > 0) {
            if (!check(fabs(energy-isolated_energy) < 1e-12, "C energy at/beyond cutoff")) return 0;
            for (int j=0; j<6; ++j)
                if (!check(fabs(forces[j]) < 1e-12, "C zero force at/beyond cutoff")) return 0;
            for (int j=0; j<9; ++j)
                if (!check(fabs(virial[j]) < 1e-12, "C zero virial at/beyond cutoff")) return 0;
        }
    }
    return 1;
}

int main(int argc, char **argv) {
    double center[3] = {0.0, 0.0, 0.0};
    double neighbors[6] = {0.3, 0.1, 0.2, -0.1, 0.5, 0.4};
    int types[2] = {1, 1}, indices[2] = {2, 3}, stat;
    double energy, old_energy, forces[9] = {0}, old_forces[9] = {0}, virial[9] = {0};
    const double h = 1e-6;
    if (argc != 2) return 2;
    accelnet_init_n2p2(argv[1], &stat);
    if (!check(stat == ACCELNET_OK, "C init")) return 1;
    accelnet_atomic_energy_and_forces(center, 1, 1, 2, neighbors, types, indices,
                                     3, &old_energy, old_forces, &stat);
    if (!check(stat == ACCELNET_OK, "C old forces")) return 1;
    accelnet_atomic_energy_and_forces_virial(center, 1, 1, 2, neighbors, types, indices,
                                            3, &energy, forces, virial, &stat);
    if (!check(stat == ACCELNET_OK && fabs(energy-old_energy) < 1e-12, "C energy compatibility")) return 1;
    for (int k=0; k<9; ++k)
        if (!check(fabs(forces[k]-old_forces[k]) < 1e-12, "C force compatibility")) return 1;
    for (int a=0; a<3; ++a) {
        for (int b=0; b<3; ++b) {
            double shifted[6], ep, em;
            memcpy(shifted, neighbors, sizeof shifted);
            for (int j=0; j<2; ++j) shifted[3*j+b] += h*neighbors[3*j+a];
            accelnet_atomic_energy(center, 1, 2, shifted, types, &ep, &stat);
            if (!check(stat == ACCELNET_OK, "C strained energy +")) return 1;
            for (int j=0; j<2; ++j) shifted[3*j+b] -= 2*h*neighbors[3*j+a];
            accelnet_atomic_energy(center, 1, 2, shifted, types, &em, &stat);
            if (!check(stat == ACCELNET_OK, "C strained energy -")) return 1;
            if (!check(fabs(virial[a+3*b]+(ep-em)/(2*h)) < 1e-7, "C virial layout/sign/strain")) return 1;
        }
    }
    {
        double saved[9];
        memcpy(saved, virial, sizeof saved);
        accelnet_atomic_energy_and_forces_virial(center, 1, 1, 2, neighbors, types, indices,
                                                3, &energy, forces, virial, &stat);
        if (!check(stat == ACCELNET_OK, "C accumulation status")) return 1;
        for (int k=0; k<9; ++k)
            if (!check(fabs(virial[k]-2*saved[k]) < 1e-12 &&
                       fabs(forces[k]-2*old_forces[k]) < 1e-12, "C additive outputs")) return 1;
    }
    if (!check_invalid_calls(1) || !check_cutoff()) return 1;
    accelnet_final(&stat);
    if (!check(stat == ACCELNET_OK, "C final before reload")) return 1;
    if (!check_invalid_calls(0)) return 1;
    /* Reinitialization must not retain an old model or its accumulator state. */
    {
        char *species[] = {"H"};
        accelnet_init(1,species,&stat);
        if (!check(stat == ACCELNET_OK, "C reinit")) return 1;
        if (!check_invalid_calls(0)) return 1;
        accelnet_load_n2p2(argv[1],&stat);
        if (!check(stat == ACCELNET_OK, "C reload")) return 1;
        if (!check_cutoff()) return 1;
    }
    accelnet_final(&stat);
    return check(stat == ACCELNET_OK, "C final") ? 0 : 1;
}
