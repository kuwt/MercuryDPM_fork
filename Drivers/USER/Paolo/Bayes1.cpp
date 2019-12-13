//
// Created by paolo on 1-10-19.
//

#include <iostream>
#include "BaseDoublePorosityStressStrain.h"
#include "DoublePorosityStressStrainA.h"

int main(){


    std::vector<Mdouble> row;
    std::vector<std::vector<Mdouble>> data;


    std::ifstream infile("smcTableA.txt");
    row.reserve(5);

    Mdouble a, b, c, d, e, f;

    while ( infile  >> a >> b >> c >> d >> e >> f)
    {
        row = {b,c,d,e,f};
        data.push_back(row);
    }

    for (int i = 0; i < data.size(); ++i) {
        std::cout << data[i][0] << "   " << data[i][1] << "   " << data[i][2] << "   " << data[i][3] << "   " << data[i][4] << std::endl;
    }

    int n = 0;

    Mdouble particleDispersity = 1.2;
    Mdouble penetrationDepthMax = 0.0179414;
    Mdouble unloadingStiffnessMaxIntra = 2.9531;
    Mdouble relativeCohesionStiffnessIntra = 1.12817;
    Mdouble slidingFriction = 0.802547;


    std::ostringstream nome1;
    nome1 << "BaseA" << n;

    BaseDoublePorosityStressStrain problem(particleDispersity,pow(10,unloadingStiffnessMaxIntra),2.0/7.0,0,1.1e-3,nome1.str());
    problem.solve();

    Mdouble boxSize = problem.getXMax();
    std::ostringstream nome2;
    nome2 << "DPBayesAValoriRiscalati_" << n << "_" << penetrationDepthMax << "_" << unloadingStiffnessMaxIntra << "_" << relativeCohesionStiffnessIntra << "_" << slidingFriction;

    DoublePorosityStressStrainA problem2(problem.clusterRadii,problem.clusterPositions,
                                         penetrationDepthMax, unloadingStiffnessMaxIntra, relativeCohesionStiffnessIntra,
                                         slidingFriction, boxSize, nome2.str());
    problem2.solve();

}