//
// Created by paolo on 12-11-19.
//
//Ciao!

#include <iostream>
#include <Mercury3D.h>

int main(){

    std::ofstream smcTable;
    std::ostringstream fileName;
    fileName << "smcTableD.txt";

    smcTable.open(fileName.str(), std::ios::out);
    smcTable << "foo\tindex\t\tExprelKcIntra\tmu\t\trelKt" << std::endl;

    RNG rand;
    rand.randomise();

    for (int i = 0; i < 72; ++i) {



        Mdouble rappCompress = rand.getRandomNumber(5, 20);//-1 + 2 * fmod(i,23.0)/23 + 2 * floor(i/23.0)/72;  //0.33334 + (100-1/3)*fmod(i,17)/17 + (100-1/3)*floor(i/17)/142;
        Mdouble relativeExpKcIntra = rand.getRandomNumber(-2, 2);//-1 + 2 * fmod(i,23.0)/23 + 2 * floor(i/23.0)/72;  //0.33334 + (100-1/3)*fmod(i,17)/17 + (100-1/3)*floor(i/17)/142;
        Mdouble mu = rand.getRandomNumber(0,0.9);//0.9*fmod(i,33)/33+0.9*floor(i/33.0)/72;
        Mdouble relTanStiff = rand.getRandomNumber(0.0, 1.0);//0.02 + 0.28*fmod(i,3.0)/3+0.28*floor(i/3.0)/72;

        smcTable << std::scientific << std::setprecision(5)
                 << 4 << "\t"
                 << i << "\t"
                 << rappCompress << "\t\t\t"
                 << relativeExpKcIntra << "\t\t\t"
                 << mu << "\t\t\t"
                 << relTanStiff
                 << std::endl;

    }

}
