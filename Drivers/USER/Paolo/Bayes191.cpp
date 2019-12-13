//
// Created by paolo on 1-10-19.
//

#include <iostream>
#include "BaseDoublePorosityStressStrainD.h"
#include "DoublePorosityStressStrainD.h"

int main(){


    std::vector<Mdouble> row;
    std::vector<std::vector<Mdouble>> data;


    std::ifstream infile("smcTableD.txt");
    row.reserve(5);

    Mdouble a, b, c, d, e;

    while ( infile  >> a >> b >> c >> d >> e)
    {
        row = {b,c,d,e};
        data.push_back(row);
    }

    for (int i = 0; i < data.size(); ++i) {
        std::cout << data[i][0] << "   " << data[i][1] << "   " << data[i][2] << "   " << data[i][3] << std::endl;
    }

    int n = 46;

    Mdouble particleDispersity = 1.2;
    Mdouble relExpKc = data[n][1];
    Mdouble slidingFriction = data[n][2];
    Mdouble relKt = data[n][3];


    std::cout << relExpKc << "   " << slidingFriction << "   " << relKt << "   " << std::endl;

    std::ostringstream nome1;
    nome1 << "BaseD" << n;

    BaseDoublePorosityStressStrainD problem(particleDispersity,pow(10,4),2.0/7.0,0,1.1e-3,nome1.str());
    problem.solve();

    Mdouble boxSize = problem.getXMax();
    std::ostringstream nome2;
    nome2 << "DPBayesD_" << n << "_" << relExpKc << "_" << slidingFriction << "_" << relKt;

    DoublePorosityStressStrainD problem2(problem.clusterRadii,problem.clusterPositions,
                                         relExpKc, slidingFriction, relKt,
                                         boxSize, boxSize, boxSize, nome2.str());
    problem2.solve();

}