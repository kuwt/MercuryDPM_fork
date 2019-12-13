//
// Created by paolo on 1-10-19.
//

#include <iostream>
#include "BaseDoublePorosityStressStrainD.h"
#include "DoublePorosityStressStrainC.h"

int main(){



    int n = 50;

    Mdouble particleDispersity = 1.2;
    Mdouble rappCompr = 0.208;
    Mdouble relExpKc = -1.73056e+00;
    Mdouble slidingFriction = 8.23733e-01;
    Mdouble relKt = 4.86594e-01;


    std::cout << rappCompr << "   "  << relExpKc << "   " << slidingFriction << "   " << relKt << "   " << std::endl;


    std::ostringstream nome1;
    nome1 << "BaseD" << n;

    BaseDoublePorosityStressStrainD problem(particleDispersity,pow(10,4),2.0/7.0,0,1.1e-3,nome1.str());
    problem.solve();

    Mdouble boxSize = problem.getXMax();
    std::ostringstream nome2;
    nome2 << "VALIDATION1_" << n << "_" << rappCompr << "_" << relExpKc << "_" << slidingFriction << "_" << relKt;

    DoublePorosityStressStrainC problem2(problem.clusterRadii,problem.clusterPositions,
                                         rappCompr, relExpKc, slidingFriction, relKt,
                                         boxSize, boxSize, boxSize, nome2.str());
    problem2.solve();

}