//
// Created by paolo on 15-9-19.
//
//Ciao!

#include <iostream>
#include "BayesianCompression.h"

int main(){
    BayesianCompression problem(1.2,1e3,0,0.0,2.4e-3, 5000,"TestVelocità_");
    problem.solve();
}
