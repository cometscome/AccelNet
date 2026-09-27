// Independent n2p2 library oracle and single-core fixed-neighbor/full timings.
#include "Prediction.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv) {
    try {
        if(argc!=5) throw std::runtime_error("usage: n2p2-reference MODEL INPUT_DATA SECONDS fixed|full");
        const double duration=std::stod(argv[3]);
        if(duration<0) throw std::runtime_error("negative duration");
        const std::string scope(argv[4]);
        if(scope!="fixed" && scope!="full") throw std::runtime_error("invalid scope");
        nnp::Prediction p; p.log.writeToStdout=false;
        std::string dir(argv[1]);
        p.fileNameSettings=dir+"/input.nn"; p.fileNameScaling=dir+"/scaling.data";
        p.formatWeightsFilesShort=dir+"/weights.%03zu.data"; p.setup();
        p.readStructureFromFile(argv[2]); auto& s=p.structure;
        p.evaluateNNP(s,true,true);
        auto eval=[&]() {
            if(scope=="full") { s.clearNeighborList(); p.evaluateNNP(s,true,true); }
            else {
                s.freeAtoms(true,p.getMaxCutoffRadius());
                p.calculateSymmetryFunctionGroups(s,true);
                p.calculateAtomicNeuralNetworks(s,true);
                p.calculateEnergy(s); p.calculateForces(s);
            }
        };
        std::cout<<std::scientific<<std::setprecision(17);
        for(int sample=0;sample<(duration>0 ? 5 : 0);sample++) {
            const auto start=std::chrono::steady_clock::now(); double elapsed=0;size_t count=0;
            do {eval();count++;elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
            while(elapsed<duration);
            std::cout<<"TIMING "<<elapsed/count<<" "<<count<<"\n";
        }
        p.convertToPhysicalUnits(s);p.addEnergyOffset(s,false);
        std::cout<<"ENERGY "<<s.energy<<"\n";
        for(const auto& a:s.atoms) std::cout<<"FORCE "<<a.f[0]<<" "<<a.f[1]<<" "<<a.f[2]<<"\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
