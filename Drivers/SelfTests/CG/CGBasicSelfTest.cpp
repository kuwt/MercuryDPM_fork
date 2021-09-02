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
//#include "StatisticsVector.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "CG/CG.h"
#include "CG/TimeSmoothedCG.h"
#include "CG/TimeAveragedCG.h"

/*! 
 * \brief Tests if the different CG templates work correctly.
 * \details The problem consists of a stack of three equal-sized particles; the 
 * bottom particle is fixed, the two particles on top are pushed towards the fixed 
 * particles by gravity. The particles are stiff, so the overlap is less than 
 * 0.03% of the particles' diameter. The problem is non-dimensionalised such that
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
 * Test that the system goes to steady state after 60% of the simulation time:
 * > set logscale y; p 'CGBasicSelfTest.ene' u 1:3; unset logscale y
 * Plot temporal data vs time-averaged data:
 * > p 'CGBasicSelfTest.LucyO.TS.stat' u 28 w l,'CGBasicSelfTest.LucyO.stat' u 28 w l
 * Plot z-data for different CG functions:
 * > set ylabel 'Density'; p 'CGBasicSelfTest.LucyZ.stat' u 5:9 ev ::2 w l, 'CGBasicSelfTest.LinearZ.stat' u 5:9 ev ::2 w l, 'CGBasicSelfTest.HeavisideZ.stat' u 5:9 ev ::2 w l, 'CGBasicSelfTest.GaussZ.stat' u 7:11 ev ::2 w l; unset ylabel
 * > set ylabel 'Stress'; p 'CGBasicSelfTest.LucyZ.stat' u 5:30 w l, 'CGBasicSelfTest.LinearZ.stat' u 5:30 w l, 'CGBasicSelfTest.HeavisideZ.stat' u 5:30 w l, 'CGBasicSelfTest.GaussZ.stat' u 7:32 w l; unset ylabel
 * 
 * > p 'CGBasicSelfTest.1.stat' u 7:32 w l, '' u 7:11 w l
 * > set xlabel 'x'; set ylabel 'z'; sp [][][0:] 'CGBasicSelfTest.2.stat' u 7:9:34 w l
 */
class CGBasicSelfTest : public Mercury3D
{
public:

    CGBasicSelfTest()
    {
        setName("CGBasicSelfTest");

        //set gravity and the species properties
        setGravity(Vec3D(-0.1, 0, -1.0));
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(6.0 / constants::pi);
        species->setCollisionTimeAndRestitutionCoefficient(.01, .1, 1.);

        //accordingly, set the time step, final time, and save count
        setTimeStep(species->getCollisionTime(1.0)/50.0);
        setSaveCount(10);
    }

    void setupInitialConditions() override {
        //define the domain
        setXMax(1.1);
        setYMax(1);
        setZMax(2.1);
        setXMin(-0.1);
        setYMin(0);
        setZMin(-0.6);

        //define the three particles stacked in Z-direction
        SphericalParticle P0;
        P0.setSpecies(speciesHandler.getLastObject());
        P0.setRadius(0.5 * sqrt(1.01));
        P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        P0.setPosition(Vec3D(0.5, 0.5, 0.5));
        particleHandler.copyAndAddObject(P0);
        P0.setPosition(Vec3D(0.6, 0.5, 1.5));
        particleHandler.copyAndAddObject(P0);
        P0.setPosition(Vec3D(0.4, 0.5, -0.5));
        P0.fixParticle();
        particleHandler.copyAndAddObject(P0);

    }
};

int main(int argc, char *argv[])
{
    //declare the DPM problem and set the name
    CGBasicSelfTest problem;
    problem.setTimeMax(0.21); //a total of 5000 time steps

    //int n1 = 1000, n2 = 100, n3=20;
    size_t n1 = 1000, n2 = 33, n3=10;
    //define different coarse-graining objects (which should all result in the same mean values)
    auto cg0 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
    cg0->setTimeMin(0.21);
    auto cg1 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::Z,CGFunctions::Gauss>());
    cg1->setN(n1);
    cg1->setWidth(0.5/3.01);
    cg1->setTimeMin(0.21);
    auto cg2 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::YZ,CGFunctions::Gauss>());
    cg2->setN(n2);
    cg2->setWidth(0.5/3.01);
    cg2->setTimeMin(0.21);
    auto cg3 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::XYZ,CGFunctions::Gauss>());
    cg3->setN(n3);
    cg3->setWidth(0.5/3.01);
    cg3->setTimeMin(0.21);
    auto cg4 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::Z,CGFunctions::Heaviside>());
    cg4->setNZ(n1);
    cg4->setWidth(0.5/2.01);
    cg4->setTimeMin(0.21);
    auto cg5 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::YZ,CGFunctions::Heaviside>());
    cg5->setN(n2);
    cg5->setWidth(0.5/2.01);
    cg5->setTimeMin(0.21);
    auto cg6 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::XYZ,CGFunctions::Heaviside>());
    cg6->setN(n3);
    cg6->setWidth(0.5/2.01);
    cg6->setTimeMin(0.21);
    auto cg7 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::Z,CGFunctions::Lucy>());
    cg7->setNZ(n1);
    cg7->setWidth(0.5/2.01);
    cg7->setTimeMin(0.21);
    auto cg8 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::YZ,CGFunctions::Lucy>());
    cg8->setN(n2);
    cg8->setWidth(0.5/2.01);
    cg8->setTimeMin(0.21);
    auto cg9 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::XYZ,CGFunctions::Lucy>());
    cg9->setN(n3);
    cg9->setWidth(0.5/2.01);
    cg9->setTimeMin(0.21);
    auto cg10 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::X,CGFunctions::Lucy>());
    cg10->setN(n1);
    cg10->setWidth(0.5/2.01);
    cg10->setTimeMin(0.21);
    auto cg11 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::X,CGFunctions::Gauss>());
    cg11->setN(n1);
    cg11->setWidth(0.5/3.01);
    cg11->setTimeMin(0.21);
    auto cg12 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::Y,CGFunctions::Lucy>());
    cg12->setN(n1);
    cg12->setWidth(0.5/2.01);
    cg12->setTimeMin(0.21);
    auto cg13 = problem.cgHandler.copyAndAddObject(CG<CGCoordinates::Y,CGFunctions::Gauss>());
    cg13->setN(n1);
    cg13->setWidth(0.5/3.01);
    cg13->setTimeMin(0.21);

    //run the simulation
    problem.solve(argc, argv);

    helpers::writeToFile("CGBasicSelfTest.1D.gnu","p 'CGBasicSelfTest.1.stat' u 2:22 w d, 'CGBasicSelfTest.4.stat' u 2:22 w d, 'CGBasicSelfTest.7.stat' u 2:22 w d");
    helpers::writeToFile("CGBasicSelfTest.2D.gnu","sp 'CGBasicSelfTest.8.stat' u 2:3:23 w d");
}

