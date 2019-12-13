//
// Created by paolo on 5-10-19.
//
//Ciao!

#include <iostream>
#include <Mercury3D.h>

int main(){

    std::ofstream smcTable;
    std::ostringstream fileName;
    fileName << "smcTableA.txt";

    smcTable.open(fileName.str(), std::ios::out);
    smcTable << "foo\tindex\t\tpenetrationDepthMax\tkeIntra\t\tExprelKcIntra\tmu" << std::endl;

    for (int i = 0; i < 35; ++i) {

        Mdouble expKeIntra = 4;//3 + 2 * fmod(i,9)/9 + 2 * floor(i/9.0)/142;
        Mdouble relativeExpKcIntra = -2 + 4 * fmod(i,5)/5 + 4 * floor(i/5.0)/35;  //0.33334 + (100-1/3)*fmod(i,17)/17 + (100-1/3)*floor(i/17)/142;
        Mdouble mu = 0.9*fmod(i,7)/7+0.9*floor(i/7.0)/35;
        Mdouble penetrationDepthMax = 0.1;//0.02 + 0.28*fmod(i,13)/13+0.28*floor(i/13.0)/142;

        smcTable << std::scientific << std::setprecision(5)
                 << 4 << "\t"
                 << i << "\t"
                 << penetrationDepthMax << "\t"
                 << expKeIntra << "\t"
                 << relativeExpKcIntra << "\t"
                 << mu
                 << std::endl;

    }

}
