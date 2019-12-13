//
// Created by paolo on 5-10-19.
//
//Ciao!

#include <iostream>
#include <Mercury3D.h>

int main(){

    std::ofstream smcTable;
    std::ostringstream fileName;
    fileName << "smcTable.txt";

    smcTable.open(fileName.str(), std::ios::out);
    smcTable << "foo\tindex\t\tpD\trelKnIntra\tkeIntra\t\tExprelKcIntra\tknInter\t\trelKt\t\tmu" << std::endl;

    for (int i = 0; i < 142; ++i) {

        Mdouble pD = 1.2 + 1.3*fmod(i,19)/19 + 1.3*floor(i/19.0)/142;
        Mdouble relativeKnIntra = 0.1 + 0.4*fmod(i,11)/11 + 0.4*floor(i/11.0)/142;
        Mdouble expKeIntra = 3 + 2 * fmod(i,15)/15 + 2 * floor(i/15.0)/142;
        Mdouble relativeExpKcIntra = -2 + 4 * fmod(i,17)/17 + 4 * floor(i/17.0)/142;  //0.33334 + (100-1/3)*fmod(i,17)/17 + (100-1/3)*floor(i/17)/142;
        Mdouble relativeKnInter = 1+99*fmod(i,13)/13+99*floor(i/13.0)/142;
        Mdouble relativeKt = 2*fmod(i,9)/9+2*floor(i/9.0)/142;
        Mdouble mu = 0.9*fmod(i,7)/7+0.9*floor(i/7.0)/142;

        smcTable << std::scientific << std::setprecision(5)
                 << 4 << "\t"
                 << i << "\t"
                 << pD << "\t"
                 << relativeKnIntra << "\t"
                 << expKeIntra << "\t"
                 << relativeExpKcIntra << "\t"
                 << relativeKnInter << "\t"
                 << relativeKt << "\t"
                 << mu
                 << std::endl;

    }

}
