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

#include "Mercury3D.h"
#include "StatisticsVector.h"
#include "Walls/InfiniteWall.h"
#include "CG/CG.h"
#include "CG/Functions/Gauss.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <Species/LinearViscoelasticSpecies.h>

/// In this file a cubic packing of 5^3 particles in a tri-axial box is created and allowed to settle under small gravity. After that Z statistics are calculated.

class CGHandlerSelfTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        double Radius = .5;

        setXMax(5);
        setYMax(5);
        setZMax(5);
        setXMin(0);
        setYMin(0);
        setZMin(0);

        SphericalParticle P0;
        P0.setSpecies(speciesHandler.getLastObject());
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                for (int k = 0; k < N; k++)
                {
                    P0.setRadius(Radius);
                    P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
                    P0.setPosition(Radius * Vec3D(1. + 2. * i, 1. + 2. * j, 1. + 2. * k));
                    particleHandler.copyAndAddObject(P0);
                }

        //set walls
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getLastObject());
        w0.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(+1, 0, 0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, -1, 0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, +1, 0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 0, +1), Vec3D(0, 0, getZMax()));
        wallHandler.copyAndAddObject(w0);
    }

    int N;

};

int main(int argc, char *argv[])
{
    CGHandlerSelfTest problem;
    problem.setName("CGHandlerSelfTest");
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());

    //set the number of particles
    problem.N = 5;
    problem.setSystemDimensions(3);
    problem.speciesHandler.getObject(0)->setDensity(6.0 / constants::pi);
    problem.setGravity(Vec3D(0., 0., -1.));
    species->setCollisionTimeAndRestitutionCoefficient(.01, .1, 1.);
    problem.setTimeStep(.0002);
    problem.setTimeMax(1.5 - 1e-10);
    problem.setSaveCount(7500);

//    CG<CGCoordinates::XY,CGFunctions::Gauss> cg0;
//    cg0.setWidth(0);
//    problem.cgHandler.copyAndAddObject(cg0);
//
//    CG<CGCoordinates::Z,CGFunctions::Gauss> cg1;
//    cg1.setNZ(100);
//    cg1.setWidth(0.1);
//    problem.cgHandler.copyAndAddObject(cg1);

    CG<CGCoordinates::XZ,CGFunctions::Gauss> cg2;
    cg2.setNX(20);
    cg2.setNZ(30);
    cg2.setWidth(0.5);
    cg2.statFile.setSaveCount(20000);
    problem.cgHandler.copyAndAddObject(cg2);

//    CG<CGCoordinates::XYZ,CGFunctions::Gauss> cg3;
//    cg3.setNX(10);
//    cg3.setNY(10);
//    cg3.setNZ(10);
//    cg3.setWidth(0.1);
//    problem.cgHandler.copyAndAddObject(cg3);

    problem.solve(argc, argv);
    //p 'CGHandlerSelfTest.1.stat' u 7:32 w l, '' u 7:11 w l
    //sp [][][0:] 'CGHandlerSelfTest.2.stat' u 7:9:34, '' u 7:9:13
}

