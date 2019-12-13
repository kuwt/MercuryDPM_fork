//
// Created by paolo on 5-10-19.
//
//Ciao!

#include <iostream>
#include <Mercury3D.h>

int main(){

    std::ofstream smcTable;
    std::ostringstream fileName;
    fileName << "smcTableB.txt";

    smcTable.open(fileName.str(), std::ios::out);
    smcTable << "foo\tindex\trelDeltaStar\t\tkeIntra\t\tExprelKcIntra\t\tmu" << std::endl;

    for (int i = 0; i < 30; ++i) {

        Mdouble relDeltaStar = 0.02 + 0.28*fmod(i,9)/9+0.28*floor(i/9.0)/30;
        Mdouble expKeIntra = 2.9531;
        Mdouble relativeExpKcIntra = 1.12817;  //0.33334 + (100-1/3)*fmod(i,17)/17 + (100-1/3)*floor(i/17)/142;
        Mdouble mu = 0.802547;

        smcTable << std::scientific << std::setprecision(5)
                 << 4 << "\t"
                 << i << "\t"
                 << relDeltaStar << "\t"
                 << expKeIntra << "\t"
                 << relativeExpKcIntra << "\t"
                 << mu
                 << std::endl;

    }


}
