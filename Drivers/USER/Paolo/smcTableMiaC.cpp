//
// Created by paolo on 21-10-19.
//
//Ciao!

#include <iostream>
#include <Mercury3D.h>

int main(){

    std::ofstream smcTable;
    std::ostringstream fileName;
    fileName << "smcTableC.txt";

    smcTable.open(fileName.str(), std::ios::out);
    smcTable << "foo\tindex\trelDeltaStar\t\tkeIntra\t\tExprelKcIntra\t\tmu" << std::endl;

    for (int i = 0; i < 72; ++i) {

        Mdouble expKeIntra = 0 + 3 * fmod(i,11)/11 + 3 * floor(i/11.0)/142;
        Mdouble relativeExpKcIntra = -2 + 4 * fmod(i,13)/13 + 4 * floor(i/13.0)/142;  //0.33334 + (100-1/3)*fmod(i,17)/17 + (100-1/3)*floor(i/17)/142;
        Mdouble mu = 0.9*fmod(i,7)/7+0.9*floor(i/7.0)/142;
        Mdouble relDeltaStar = 0.1231 + 0.1*fmod(i,9)/9+0.1*floor(i/9.0)/142;

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
