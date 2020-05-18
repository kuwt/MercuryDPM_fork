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
#include <cmath>
#include <iostream>
#include <iomanip>
#include <Species/LinearViscoelasticSpecies.h>
#include <CG/Functions/Heaviside.h>
#include <Species/LinearViscoelasticReversibleAdhesiveSpecies.h>
#include <CG/TimeAveragedCG.h>
#include "MercuryTime.h"
#include "CG/CG.h"
using constants::pi;
using mathsFunc::cubic;

/// In this file a cubic packing of 5^3 particles in a tri-axial box is created and allowed to settle under small gravity. After that Z statistics are calculated.

struct MercuryCGSelfTest : public Mercury3D
{
    explicit MercuryCGSelfTest(const unsigned n=5)
    {
        setMin({0,0,0});
        setMax(Vec3D(n,1,1));
        setGravity({0,0,0});
        setName("MercuryCGSelfTest");
        setTimeStep(1);
        setTimeMax(2);
        setSaveCount(1);

        LinearViscoelasticReversibleAdhesiveSpecies s;
        s.setDensity(6./pi); //such that density=1
        s.setCollisionTimeAndRestitutionCoefficient(1, 1, 1);
        s.setAdhesionStiffness(s.getStiffness());
        s.setAdhesionForceMax(1.0);//such that force = 1 at distance 1 (contactStress=-1)
        speciesHandler.copyAndAddObject(s);
        speciesHandler.copyAndAddObject(s);

        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(0.5);
        for (Mdouble x = getXMin()+0.5; x < getXMax(); x++)
            for (Mdouble y = getYMin()+0.5; y < getYMax(); y++)
                for (Mdouble z = getZMin()+0.5; z < getZMax(); z++)
                {
                    p.setPosition({x,y,z});
                    particleHandler.copyAndAddObject(p);
                }
        particleHandler.getLastObject()->setSpecies(speciesHandler.getLastObject());

        //set walls
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1, 0, 0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
//        w0.set(Vec3D(0, -1, 0), Vec3D(0, getYMin(), 0));
//        wallHandler.copyAndAddObject(w0);
//        w0.set(Vec3D(0, 1, 0), Vec3D(0, getYMax(), 0));
//        wallHandler.copyAndAddObject(w0);
//        w0.set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
//        wallHandler.copyAndAddObject(w0);
//        w0.set(Vec3D(0, 0, 1), Vec3D(0, 0, getZMax()));
//        wallHandler.copyAndAddObject(w0);
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    unsigned n = 2;
    MercuryCGSelfTest problem(n);
    problem.solve();

    Mercury3D cg;
    auto c0 = cg.cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::XYZ,CGFunctions::Heaviside>());
    //c0->setX(-0.5,n+0.5);
    c0->setNX(250);
    c0->setNY(1);
    c0->setNZ(1);
    c0->setWidth(0.3);
    //c0->selectSpecies(0);

    //cg.cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
    cg.cgHandler.restartAndEvaluateDataFiles("MercuryCGSelfTest");

    helpers::more(cg.cgHandler.getLastObject()->statFile.getName(),3);

    logger(INFO,"Run 'sp 'MercuryCGSelfTest.0.stat' u ($4==0.5?$2:NaN):3:16' in gnuplot to view the output");

    helpers::writeToFile("MercuryCGSelfTest.gnu","set multiplot layout 1,2\n"
     "set xlabel 'x'\n"
     //"set ylabel 'y'\n"
     "set ylabel rotate\n"
     "unset key\n"
     "set ylabel 'Density at z=0.5'\n"
     //"sp 'MercuryCGSelfTest.0.stat' u ($4==0.5?$2:NaN):3:6\n"
     "p 'MercuryCGSelfTest.0.stat' u ($4==0.5&&$3==0.5?$2:NaN):6 w lp\n"
     "unset key\n"
     "set ylabel 'Normal Contact Stress in x at z=0.5'\n"
     //"sp 'MercuryCGSelfTest.0.stat' u ($4==0.5?$2:NaN):3:16\n"
     "p 'MercuryCGSelfTest.0.stat' u ($4==0.5&&$3==0.5?$2:NaN):16 w lp\n"
     "unset multiplot\n");
}

