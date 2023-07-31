//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include "Helpers/Helpers.h"

typedef Eigen::SparseMatrix<double> mat;
typedef Eigen::Triplet<double> T;
std::vector<T> tripletList;

//solved the heat equation dT/dt=1*(d^2T)/(dx^2)
int main(int argc, char** argv){
    double L = 1;
    int varnum = 5;
    double dx = L/(varnum-1);
    double dt = pow(dx,2);//does not need to be this, but keeps it anyway
    double A = dt/2.0/pow(dx,2);
    double B = -(1+A*2);


    //construct coefficient sparse matrix
    //diagonal elements
    for(int i=0; i<varnum; i++) {
        tripletList.push_back(T(i, i, B));
    }

    //right-diagonal elements
    for(int j=1; j<varnum; j++){
        tripletList.push_back(T(j-1,j,A));
    }
    //left-diagonal elements
    for(int k=0; k<varnum-1; k++){
        tripletList.push_back(T(k+1, k, A));
    }

    //distribute the elements to the matrix
    mat C(varnum,varnum);
    C.setFromTriplets(tripletList.begin(), tripletList.end());

    //construct RHS vector
    Eigen::VectorXd T(varnum), K(varnum);
    for(int i=0; i<varnum; i++)
        T(i) = 0;

    for(int j=1; j<varnum-1; j++){
        K(j) = -T(j)-A*(T(j+1)-2*T(j)+T(j-1));
    }
    K(0) = -T(0)-A*(T(1)-2*T(0));
    K(varnum-1) = -T(varnum-1)-A*(1-2*T(varnum-1)+T(varnum-2))-A;


    //iterate with time
    Eigen::ConjugateGradient<mat> solver;
    for(double t=0; t<dt*100; t+=dt){
        solver.compute(C);
        T = solver.solve(K);

        for(int j=1; j<varnum-1; j++){
            K(j) = -T(j)-A*(T(j+1)-2*T(j)+T(j-1));
        }
        K(0) = -T(0)-A*(T(1)-2*T(0));
        K(varnum-1) = -T(varnum-1)-A*(1-2*T(varnum-1)+T(varnum-2))-A;

        //output the values at the grid points to the console
        std::cout << "t=" << t << "s" << "\t" << T(0) << "\t" << T((varnum-1)*1/4) << "\t" << T((varnum-1)/2) << "\t" << T((varnum-1)*3/4) << "\t" << T(varnum-1)<< std::endl;
        
        
        //check last timestep
        if (6.1875 - t  <= 1e-5)
        {
            helpers::check(T(0),0.166667,1e-5, "Gridpoint 0");
            helpers::check(T((varnum-1)*1/4),0.333333,1e-5, "Gridpoint 1/4");
            helpers::check(T((varnum-1)*2/4),0.5,1e-5, "Gridpoint 2/4");
            helpers::check(T((varnum-1)*3/4),0.666667,1e-5, "Gridpoint 3/4");
            helpers::check(T(varnum-1),0.833333,1e-5, "Gridpoint 1");
        }
    }

}





