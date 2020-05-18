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
#include "CG/TimeSmoothedCG.h"
#include "CG/TimeAveragedCG.h"
#include "CG/Functions/Lucy.h"
#include "CG/Functions/Linear.h"
#include "CG/Functions/Heaviside.h"
#include "CG/Functions/Gauss.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <Species/LinearViscoelasticSpecies.h>

/*! 
 * \brief Tests if the different CG templates work correctly.
 * \details The problem consists of a stack of three equal-sized particles; the 
 * bottom particle is fixed, the two particles on top are pushed towards the fixed 
 * particles by gravity. The particles are stiff, so the overlap is less than 
 * 0.03% of the particles' diameter. The proplem is non-dimensionalised such that 
 * particle diameter, mass, and gravity are unity. 
 * 
 * This makes the bulk statistics very easy to compute: The bulk mass is 2 
 * (since the fixed particle is treated as external/boundary and does not 
 * contribute to the bulk, just like a wall), the bulk stress is almost 2 (force 
 * times length of branch vector almost equals unity for each of the two 
 * contacts, except for a small correction due to the overlap), and the 
 * interaction force with the boundary (the force acting on the fixed particle) 
 * is 2. Since the domain size is 1x1x2.5, the averages for density, stress and 
 * interaction force density are all 0.8 (to be exact, stress is 0.799984).
 * 
 * Several CG options are tested:
 * - computing time-averaged statistics for the last half of the simulation, 
 *   with O and XYZ coordinates using Lucy CG functions.
 * - computing statistics for the full duration of the simulation 
 *   with O and XYZ coordinates, using Lucy CG functions.
 * - computing time-smoothed statistics for the full duration of the simulation 
 *   with O and XYZ coordinates using Lucy CG functions.
 * - computing statistics for the last time step with Z, XZ, XYZ coordinates, 
 *   using Gaussian, Heaviside, Linear, and Lucy CG functions.
 * The CG widths are chosen such that the coarse-grained fields are fully in the 
 * domain.
 * 
 * The output can be tested with the following gnuplot commands:
 * 
 * > set xlabel 'r/x'; p 'CGHandlerRadialSelfTest.GaussX.stat' u 7:24 w lp, 'CGHandlerRadialSelfTest.GaussR.stat' u 7:($7*$24) w l
 */
class CGHandlerSelfTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        setName("CGHandlerRadialSelfTest");

        //species properties
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(6.0 / constants::pi);
        species->setCollisionTimeAndRestitutionCoefficient(.01, .1, 1);
        setTimeStep(.02*.01);

        //define the domain
        setMin({-1,-1,-1});
        setMax({3,3,3});

        //set gravity
        setGravity(Vec3D::getUnitVector(Vec3D(-1,-1,-0)));

        //define the three particles
        SphericalParticle p;
        p.setSpecies(speciesHandler.getLastObject());
        p.setRadius(0.5 * getGravity().getLength());
        p.setPosition(-1.5 * getGravity());
        particleHandler.copyAndAddObject(p);
        p.setPosition(-2.5 * getGravity());
        particleHandler.copyAndAddObject(p);
        p.setPosition(-0.5 * getGravity());
        p.fixParticle();
        particleHandler.copyAndAddObject(p);
    }
};

int main(int argc, char *argv[])
{
    //declare the DPM problem and set the name
    CGHandlerSelfTest problem;
    problem.setTimeMax(1.0);
    problem.setSaveCount(999999);

    problem.cgHandler.copyAndAddObject(CG<CGCoordinates::X,CGFunctions::Gauss>());
    problem.cgHandler.copyAndAddObject(CG<CGCoordinates::R,CGFunctions::Gauss>());
//    problem.cgHandler.copyAndAddObject(CG<CGCoordinates::XZ,CGFunctions::Gauss>());
//    problem.cgHandler.copyAndAddObject(CG<CGCoordinates::RZ,CGFunctions::Gauss>());

    for (auto cg : problem.cgHandler)
    {
        cg->setN(100);
        cg->setWidth(0.5/3.01);
        cg->setTimeMin(problem.getTimeMax());
    }

    //run the simulation
    problem.solve(argc, argv);

    helpers::writeToFile("CGHandlerRadialSelfTest.0.gnu","p 'CGHandlerRadialSelfTest.0.stat' u 2:($3*16), 'CGHandlerRadialSelfTest.1.stat' u 2:($3*2*pi*$2*4)");
    helpers::writeToFile("CGHandlerRadialSelfTest.1.gnu","p 'CGHandlerRadialSelfTest.0.stat' u 2:($14*16), 'CGHandlerRadialSelfTest.1.stat' u 2:($14*2*pi*$2*4)");

    //helpers::gnuplot("p 'CGHandlerRadialSelfTest.0.stat' u 2:($3*3), 'CGHandlerRadialSelfTest.1.stat' u 2:($3*2*pi*$2)\n");
}

