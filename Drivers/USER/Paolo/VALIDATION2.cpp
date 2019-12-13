//
// Created by paolo on 1-10-19.
//

#include <iostream>
#include "BaseDoublePorosityStressStrainD.h"
#include "DoublePorosityStressStrainC.h"

int main(){


    int n = 41;

    Mdouble particleDispersity = 1.2;
    Mdouble rappCompr = 0.170;
    Mdouble relExpKc = -1.78168e+00;
    Mdouble slidingFriction = 7.91489e-01;
    Mdouble relKt = 5.55923e-01;


    std::cout << rappCompr << "   "  << relExpKc << "   " << slidingFriction << "   " << relKt << "   " << std::endl;


    std::ostringstream nome1;
    nome1 << "BaseD" << n;

    BaseDoublePorosityStressStrainD problem(particleDispersity,pow(10,4),2.0/7.0,0,1.1e-3,nome1.str());
    problem.solve();

    Mdouble boxSize = problem.getXMax();
    std::ostringstream nome2;
    nome2 << "VALIDARTION2_" << n << "_" << rappCompr << "_" << relExpKc << "_" << slidingFriction << "_" << relKt;

    DoublePorosityStressStrainC problem2(problem.clusterRadii,problem.clusterPositions,
                                         rappCompr, relExpKc, slidingFriction, relKt,
                                         boxSize, boxSize, boxSize, nome2.str());
    problem2.solve();

}