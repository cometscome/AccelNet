#include "accelnet_target.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define REQUIRE(x) do { if (!(x)) { fprintf(stderr,"line %d: %s (%s)\n",__LINE__,#x,error); return 1; } } while (0)
int main(int argc, char **argv) {
    if (argc != 3) return 2;
    const char *paths[2] = {argv[1],argv[2]};
    char error[512] = {0};
    int species[2]={1,2}, centers[2]={1,2}, offsets[3]={1,2,3}, indices[2]={2,1};
    double dr[6]={2,0,0,-2,0,0}, reference[2], rf[6]={0}, rw[9]={0};
    double cutoff;
    void *first=NULL, *second=NULL;
    REQUIRE(accelnet_target_create(2,paths,0,0,&cutoff,&first,error)==0);
    REQUIRE(first && cutoff>0);
    REQUIRE(accelnet_target_compute(first,2,2,2,species,centers,offsets,indices,dr,reference,rf,rw,error)==0);
    for (int iteration=0;iteration<8;iteration++) {
        double e[2],f[6]={0},w[9]={0};
        REQUIRE(accelnet_target_create_modes(2,paths,0,iteration%3,iteration%4,&cutoff,&second,error)==0);
        REQUIRE(second && second!=first);
        offsets[0]=0;
        REQUIRE(accelnet_target_compute(second,2,2,2,species,centers,offsets,indices,dr,e,f,w,error)!=0);
        REQUIRE(strstr(error,"offset")!=NULL);
        offsets[0]=1;
        REQUIRE(accelnet_target_compute(second,2,2,2,species,centers,offsets,indices,dr,e,f,w,error)==0);
        for (int i=0;i<2;i++) REQUIRE(isfinite(e[i]) && fabs(e[i]-reference[i])<1e-10);
        for (int i=0;i<6;i++) REQUIRE(isfinite(f[i]) && fabs(f[i]-rf[i])<1e-10);
        for (int i=0;i<9;i++) REQUIRE(isfinite(w[i]) && fabs(w[i]-rw[i])<1e-10);
        accelnet_target_destroy(second); second=NULL;
    }
    accelnet_target_destroy(first);
    accelnet_target_destroy(NULL);
    REQUIRE(accelnet_target_create(2,paths,0,99,&cutoff,&second,error)!=0 && !second);
    REQUIRE(accelnet_target_create_modes(2,paths,0,0,99,&cutoff,&second,error)!=0 && !second);
    const char *reversed[2]={argv[2],argv[1]};
    REQUIRE(accelnet_target_create(2,reversed,0,0,&cutoff,&second,error)==0 && second);
    double e[2],f[6]={0},w[9]={0};
    REQUIRE(accelnet_target_compute(second,2,2,2,species,centers,offsets,indices,dr,e,f,w,error)==0);
    for (int i=0;i<2;i++) REQUIRE(fabs(e[i]-reference[i])<1e-10);
    for (int i=0;i<6;i++) REQUIRE(fabs(f[i]-rf[i])<1e-10);
    for (int i=0;i<9;i++) REQUIRE(fabs(w[i]-rw[i])<1e-10);
    accelnet_target_destroy(second); second=NULL;
    const char *duplicate[2]={argv[1],argv[1]};
    REQUIRE(accelnet_target_create(2,duplicate,0,0,&cutoff,&second,error)!=0 && !second);
    FILE *bad=fopen("target-invalid.nn.ascii","w");
    REQUIRE(bad!=NULL);
    fputs("4\n56\n",bad); fclose(bad);
    const char *invalid[2]={"target-invalid.nn.ascii",argv[1]};
    REQUIRE(accelnet_target_create(2,invalid,0,0,&cutoff,&second,error)!=0 && !second);
    remove("target-invalid.nn.ascii");
    puts("GPU C API: independent instances, validation, cleanup passed");
    return 0;
}
