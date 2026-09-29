#include "accelnet_target.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char error[512];
#define REQUIRE(x) do { if (!(x)) { fprintf(stderr,"line %d: %s (%s)\n",__LINE__,#x,error); exit(1); } } while (0)
static const double positions[12]={0,0,0, 1.1,.2,.1, .3,1.3,-.2, 1.5,1.1,.4};
static int centers[4]={1,2,3,4}, offsets[5]={1,4,7,10,13}, indices[12];
static double dr[36];

static void close_values(const double *actual,const double *reference,int n) {
    for (int i=0;i<n;i++) REQUIRE(isfinite(actual[i]) &&
        fabs(actual[i]-reference[i]) < 2e-10+2e-10*fabs(reference[i]));
}
static void check(void *handle,int nspecies,double cutoff,const char *reference) {
    int n,species[4];
    double rc,expected[25],e[4],f[12]={0},w[9]={0};
    char expected_symbol[17],symbol[17];
    FILE *file=fopen(reference,"r");
    REQUIRE(file && fscanf(file,"%d",&n)==1 && n==nspecies);
    for (int i=1;i<=n;i++) {
        REQUIRE(fscanf(file,"%16s",expected_symbol)==1);
        REQUIRE(accelnet_target_get_species(handle,i,sizeof(symbol),symbol,error)==0);
        REQUIRE(strcmp(symbol,expected_symbol)==0);
    }
    REQUIRE(fscanf(file,"%lf",&rc)==1 && fabs(rc-cutoff)<1e-12);
    for (int i=0;i<25;i++) REQUIRE(fscanf(file,"%lf",&expected[i])==1);
    fclose(file);
    for (int i=0;i<4;i++) species[i]=1+i%nspecies;
    REQUIRE(accelnet_target_compute(handle,4,4,12,species,centers,offsets,indices,dr,e,f,w,error)==0);
    close_values(e,expected,4); close_values(f,expected+4,12); close_values(w,expected+16,9);
    /* Reuse the resident workspace and accumulate row subsets, as in XMPI. */
    memset(f,0,sizeof(f)); memset(w,0,sizeof(w));
    for (int i=0;i<4;i++)
        REQUIRE(accelnet_target_compute(handle,4,1,12,species,centers+i,offsets+i,indices,dr,e+i,f,w,error)==0);
    close_values(e,expected,4); close_values(f,expected+4,12); close_values(w,expected+16,9);
    REQUIRE(accelnet_target_get_species(handle,0,sizeof(symbol),symbol,error)!=0 && symbol[0]=='\0');
    REQUIRE(accelnet_target_get_species(handle,n+1,sizeof(symbol),symbol,error)!=0);
    REQUIRE(accelnet_target_get_species(handle,1,1,symbol,error)!=0 && symbol[0]=='\0');
}
int main(int argc,char **argv) {
    REQUIRE(argc==4);
    int device=strcmp(argv[3],"host")==0 ? ACCELNET_TARGET_HOST : 0;
    const int versions[3]={0,1,10};
    const char *fixtures[4]={"n2p2","n2p2-per-element","n2p2-per-element-depth","n2p2-virial-angular"};
    char hpath[4096],opath[4096],reference[4096],directory[4096];
    snprintf(hpath,sizeof(hpath),"%s/H.cheb.nn.ascii",argv[1]);
    snprintf(opath,sizeof(opath),"%s/O.cheb.nn.ascii",argv[1]);
    const char *paths[2]={hpath,opath};
    const char *reversed[2]={opath,hpath};
    double cutoff;
    void *handle=NULL,*other=NULL;
    int edge=0,nspecies;
    for (int i=0;i<4;i++) for (int j=0;j<4;j++) if (i!=j) {
        indices[edge]=j+1;
        for (int axis=0;axis<3;axis++) dr[3*edge+axis]=positions[3*j+axis]-positions[3*i+axis];
        edge++;
    }
    for (int v=0;v<3;v++) for (int mode=0;mode<3;mode++) {
        REQUIRE(accelnet_target_create_versioned(2,mode%2 ? reversed : paths,device,
                                                 versions[v],mode,0,&cutoff,&handle,error)==0);
        snprintf(reference,sizeof(reference),"%s/cheb-%d.ref",argv[1],versions[v]);
        check(handle,2,cutoff,reference);
        /* A second live handle must not overwrite the first model/version. */
        double other_cutoff;
        REQUIRE(accelnet_target_create_versioned(2,paths,device,versions[(v+1)%3],mode,0,
                                                 &other_cutoff,&other,error)==0);
        REQUIRE(other!=handle);
        check(handle,2,cutoff,reference);
        accelnet_target_destroy(other); other=NULL;
        accelnet_target_destroy(handle); handle=NULL;
    }
    /* Both older constructors still mean version zero. */
    snprintf(reference,sizeof(reference),"%s/cheb-0.ref",argv[1]);
    REQUIRE(accelnet_target_create(2,paths,device,0,&cutoff,&handle,error)==0);
    check(handle,2,cutoff,reference); accelnet_target_destroy(handle); handle=NULL;
    REQUIRE(accelnet_target_create_modes(2,paths,device,1,0,&cutoff,&handle,error)==0);
    check(handle,2,cutoff,reference); accelnet_target_destroy(handle); handle=NULL;
    for (int i=0;i<4;i++) for (int mode=0;mode<4;mode++) {
        snprintf(directory,sizeof(directory),"%s/%s",argv[2],fixtures[i]);
        REQUIRE(accelnet_target_create_n2p2(directory,device,mode,&nspecies,&cutoff,&handle,error)==0);
        snprintf(reference,sizeof(reference),"%s/%s.ref",argv[1],fixtures[i]);
        check(handle,nspecies,cutoff,reference);
        accelnet_target_destroy(handle); handle=NULL;
    }
    REQUIRE(accelnet_target_create_versioned(2,paths,device,2,0,0,&cutoff,&handle,error)!=0);
    REQUIRE(!handle && cutoff==0 && strstr(error,"version"));
    REQUIRE(accelnet_target_create_versioned(2,paths,device,0,99,0,&cutoff,&handle,error)!=0 && !handle);
    REQUIRE(accelnet_target_create_versioned(2,paths,device,0,0,99,&cutoff,&handle,error)!=0 && !handle);
    const char *duplicate[2]={hpath,hpath};
    REQUIRE(accelnet_target_create_versioned(2,duplicate,device,0,0,0,&cutoff,&handle,error)!=0 && !handle);
    REQUIRE(strstr(error,"Duplicate"));
    REQUIRE(accelnet_target_create_n2p2(NULL,device,0,&nspecies,&cutoff,&handle,error)!=0);
    REQUIRE(!handle && nspecies==0 && cutoff==0 && error[0]);
    REQUIRE(accelnet_target_create_n2p2("",device,0,&nspecies,&cutoff,&handle,error)!=0 && !handle);
    snprintf(directory,sizeof(directory),"%s/missing",argv[1]);
    REQUIRE(accelnet_target_create_n2p2(directory,device,0,&nspecies,&cutoff,&handle,error)!=0 && !handle);
    REQUIRE(accelnet_target_create_n2p2(argv[2],device,4,&nspecies,&cutoff,&handle,error)!=0 && !handle);
    REQUIRE(accelnet_target_create_n2p2(argv[2],-2,0,&nspecies,&cutoff,&handle,error)!=0 && !handle);
    char symbol[17];
    REQUIRE(accelnet_target_get_species(NULL,1,sizeof(symbol),symbol,error)!=0 && symbol[0]=='\0');
    accelnet_target_destroy(NULL);
    puts("C target loading: versions 0/1/10, n2p2, E/F/virial, row subsets and lifecycle passed");
    return 0;
}
