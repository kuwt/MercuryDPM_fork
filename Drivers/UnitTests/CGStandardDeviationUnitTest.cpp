//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#include <Walls/InfiniteWall.h>
#include <CG/Functions/Linear.h>
#include "Mercury3D.h"
#include "CG/CG.h"
#include "Species/LinearViscoelasticSpecies.h"
using namespace constants;

int main()
{
    logger(INFO," Checking setting the standard deviation of cg functions");

    {
        CG<CGCoordinates::X, CGFunctions::Heaviside> X;
        X.setStandardDeviation(1);
        helpers::check(X.getFunction().getCutoff(),constants::sqrt_3,1e-11, "Cutoff of 1D Heaviside with unit variance");

        CG<CGCoordinates::XY, CGFunctions::Heaviside> XY;
        XY.setStandardDeviation(1);
        helpers::check(XY.getFunction().getCutoff(),constants::sqrt_2,1e-11, "Cutoff of 2D Heaviside with unit variance");

        CG<CGCoordinates::XYZ, CGFunctions::Heaviside> XYZ;
        XYZ.setStandardDeviation(1);
        helpers::check(XYZ.getFunction().getCutoff(),sqrt(5./3.),1e-11, "Cutoff of 3D Heaviside with unit variance");
    }

    {
        CG<CGCoordinates::X, CGFunctions::Lucy> X;
        X.setStandardDeviation(1);
        helpers::check(X.getFunction().getCutoff(),sqrt(10.5),1e-11, "Cutoff of 1D Lucy with unit variance");

        CG<CGCoordinates::XY, CGFunctions::Lucy> XY;
        XY.setStandardDeviation(1);
        helpers::check(XY.getFunction().getCutoff(),sqrt(5.6),1e-11, "Cutoff of 2D Lucy with unit variance");

        CG<CGCoordinates::XYZ, CGFunctions::Lucy> XYZ;
        XYZ.setStandardDeviation(1);
        helpers::check(XYZ.getFunction().getCutoff(),2,1e-11, "Cutoff of 3D Lucy with unit variance");
    }

    {
        CG<CGCoordinates::X,CGFunctions::Gauss> X;
        X.setStandardDeviation(1);
        helpers::check(X.getFunction().getCutoff(),3,1e-11, "Cutoff of 1D Gauss with unit variance");

        CG<CGCoordinates::XY,CGFunctions::Gauss> XY;
        XY.setStandardDeviation(1);
        helpers::check(XY.getFunction().getCutoff(),3/sqrt(2),1e-11, "Cutoff of 2D Gauss with unit variance");

        CG<CGCoordinates::XYZ,CGFunctions::Gauss> XYZ;
        XYZ.setStandardDeviation(1);
        helpers::check(XYZ.getFunction().getCutoff(),3/sqrt(3),1e-11, "Cutoff of 3D Gauss with unit variance");
    }



    {
        CG<CGCoordinates::X, CGFunctions::Heaviside> X;
        X.setRadius(1);
        helpers::check(X.getFunction().getCutoff(),constants::sqrt_3*sqrt(.2),1e-11, "Cutoff of 1D Heaviside with unit radius-equivalent");

        CG<CGCoordinates::XY, CGFunctions::Heaviside> XY;
        XY.setRadius(1);
        helpers::check(XY.getFunction().getCutoff(),constants::sqrt_2*sqrt(.4),1e-11, "Cutoff of 2D Heaviside with unit radius-equivalent");

        CG<CGCoordinates::XYZ, CGFunctions::Heaviside> XYZ;
        XYZ.setRadius(1);
        helpers::check(XYZ.getFunction().getCutoff(),sqrt(5./3.)*sqrt(.6),1e-11, "Cutoff of 3D Heaviside with unit radius-equivalent");
    }

    {
        CG<CGCoordinates::X, CGFunctions::Linear> X;
        X.setRadius(1);
        helpers::check(X.getFunction().getCutoff(),sqrt(6)*sqrt(.2),1e-11, "Cutoff of 1D Linear with unit radius-equivalent");

        CG<CGCoordinates::XY, CGFunctions::Linear> XY;
        XY.setRadius(1);
        helpers::check(XY.getFunction().getCutoff(),sqrt(10./3)*sqrt(.4),1e-11, "Cutoff of 2D Linear with unit radius-equivalent");

        CG<CGCoordinates::XYZ, CGFunctions::Linear> XYZ;
        XYZ.setRadius(1);
        helpers::check(XYZ.getFunction().getCutoff(),sqrt(5./2)*sqrt(.6),1e-11, "Cutoff of 3D Linear with unit radius-equivalent");
    }

    {
        CG<CGCoordinates::X, CGFunctions::Lucy> X;
        X.setRadius(1);
        helpers::check(X.getFunction().getCutoff(),sqrt(10.5)*sqrt(.2),1e-11, "Cutoff of 1D Lucy with unit radius-equivalent");

        CG<CGCoordinates::XY, CGFunctions::Lucy> XY;
        XY.setRadius(1);
        helpers::check(XY.getFunction().getCutoff(),sqrt(5.6)*sqrt(.4),1e-11, "Cutoff of 2D Lucy with unit radius-equivalent");

        CG<CGCoordinates::XYZ, CGFunctions::Lucy> XYZ;
        XYZ.setRadius(1);
        helpers::check(XYZ.getFunction().getCutoff(),2*sqrt(.6),1e-11, "Cutoff of 3D Lucy with unit radius-equivalent");
    }

    {
        CG<CGCoordinates::X,CGFunctions::Gauss> X;
        X.setRadius(1);
        helpers::check(X.getFunction().getCutoff(),3*sqrt(.2),1e-11, "Cutoff of 1D Gauss with unit radius-equivalent");

        CG<CGCoordinates::XY,CGFunctions::Gauss> XY;
        XY.setRadius(1);
        helpers::check(XY.getFunction().getCutoff(),3/sqrt(2)*sqrt(.4),1e-11, "Cutoff of 2D Gauss with unit radius-equivalent");

        CG<CGCoordinates::XYZ,CGFunctions::Gauss> XYZ;
        XYZ.setRadius(1);
        helpers::check(XYZ.getFunction().getCutoff(),3/sqrt(3)*sqrt(.6),1e-11, "Cutoff of 3D Gauss with unit radius-equivalent");
    }
    return 0;
}

//Matlab script to compute standard deviations:

//%% Goal compute reative standard deviation, i.e.
//% var = int r^2 phi(r,w=1) dr
//syms r n real
//assume(in(n,'integer') & n>0)
//format compact
//%% Gauss function in 1D/2D/3D
//var1  = eval(int(r^2*2*exp(-r^2/2)/sqrt(2*pi),r,0,inf)); %1D
//var2  = eval(int(r^2*(2*pi*r)*exp(-r^2/2)/(2*pi),r,0,inf)); %2D
//var3  = eval(int(r^2*(4*pi*r^2)*exp(-r^2/2)/(2*pi)^1.5,r,0,inf)); %3D
//vol1  = eval(int(2*exp(-r^2/2)/sqrt(2*pi),r,0,inf)); %1D
//vol2  = eval(int((2*pi*r)*exp(-r^2/2)/(2*pi),r,0,inf)); %2D
//vol3  = eval(int((4*pi*r^2)*exp(-r^2/2)/(2*pi)^1.5,r,0,inf)); %3D
//
//%% Polynomial x^n, x<1 in 1D/2D/3D
//var1  = eval(int(2*r^2*r^n,r,0,1)) %1D
//var2  = eval(int((2*pi*r)*r^2*r^n,r,0,1)) %2D
//var3  = eval(int((4*pi*r^2)*r^2*r^n,r,0,1)) %3D
//vol1  = eval(int(2*r^n,r,0,1)); %1D
//vol2  = eval(int((2*pi*r)*r^n,r,0,1)); %2D
//vol3  = eval(int((4*pi*r^2)*r^n,r,0,1)); %3D
//
//%% Real particle in 1D/2D/3D
//vol = 4/3*pi;
//var1  = eval(int(r^2*(2)*pi*(1-r^2)/vol,r,0,1)) %1D
//var2  = eval(int(r^2*(2*pi*r)*2*sqrt(1-r^2)/vol,r,0,1)) %2D
//var3  = eval(int(r^2*(4*pi*r^2)/vol,r,0,1)) %3D
//vol1  = eval(int((2)*pi*(1-r^2)/vol,r,0,1)); %1D
//vol2  = eval(int((2*pi*r)*2*sqrt(1-r^2)/vol,r,0,1)); %2D
//vol3  = eval(int((4*pi*r^2)/vol,r,0,1)); %3D

