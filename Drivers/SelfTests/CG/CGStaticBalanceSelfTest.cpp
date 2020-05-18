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
 * Test that the system goes to steady state after 60% of the simulation time:
 * > set logscale y; p 'CGStaticBalanceSelfTest.ene' u 1:3; unset logscale y
 * Plot temporal data vs time-averaged data:
 * > p 'CGStaticBalanceSelfTest.LucyO.TS.stat' u 28 w l,'CGStaticBalanceSelfTest.LucyO.stat' u 28 w l
 * Plot z-data for different CG functions:
 * > set ylabel 'Density'; p 'CGStaticBalanceSelfTest.LucyZ.stat' u 5:9 ev ::2 w l, 'CGStaticBalanceSelfTest.LinearZ.stat' u 5:9 ev ::2 w l, 'CGStaticBalanceSelfTest.HeavisideZ.stat' u 5:9 ev ::2 w l, 'CGStaticBalanceSelfTest.GaussZ.stat' u 7:11 ev ::2 w l; unset ylabel
 * > set ylabel 'Stress'; p 'CGStaticBalanceSelfTest.LucyZ.stat' u 5:30 w l, 'CGStaticBalanceSelfTest.LinearZ.stat' u 5:30 w l, 'CGStaticBalanceSelfTest.HeavisideZ.stat' u 5:30 w l, 'CGStaticBalanceSelfTest.GaussZ.stat' u 7:32 w l; unset ylabel
 * 
 * > p 'CGStaticBalanceSelfTest.1.stat' u 7:32 w l, '' u 7:11 w l
 * > set xlabel 'x'; set ylabel 'z'; sp [][][0:] 'CGStaticBalanceSelfTest.2.stat' u 7:9:34 w l
 */
class CGStaticBalanceSelfTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        //set gravity vector
        setGravity(Vec3D(0., 0., -1.));

        //set basic contact force model
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(6.0 / constants::pi);
        species->setCollisionTimeAndRestitutionCoefficient(.01, .1, 1.);

        //accordingly, set the time step, final time, and save count
        setTimeStep(.0002);
        setTimeMax(1.0); //a total of 5000 time steps, steady after t=0.6
        setSaveCount(20);

        //define particle size
        Mdouble d = 1.0;

        //define the domain
        setXMax(d);
        setYMax(d);
        setZMax(2 * d);
        setXMin(0);
        setYMin(0);
        setZMin(-0.5 * d);

        //define the three particles stacked in z-direction
        SphericalParticle P0;
        P0.setSpecies(speciesHandler.getLastObject());
        P0.setRadius(0.5 * d);
        P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        P0.setPosition(d * Vec3D(0.5, 0.5, 0.5));
        particleHandler.copyAndAddObject(P0);
        P0.setPosition(d * Vec3D(0.5, 0.5, 1.5));
        particleHandler.copyAndAddObject(P0);
        P0.setPosition(d * Vec3D(0.5, 0.5, -0.5));
        P0.fixParticle();
        particleHandler.copyAndAddObject(P0);

    }
};


void setCGHandler(DPMBase& problem)
{
    //define the different coarse-graining objects
    
    //case 0: create global average rho(t)
    CG<CGCoordinates::O> cg0;
    cg0.statFile.setSaveCount(20);
    cg0.statFile.setName(problem.getName() + ".LucyO.stat");
    problem.cgHandler.copyAndAddObject(cg0);

//    TimeSmoothedCG<CGCoordinates::O> cg1;
//    cg1.setWidthTime(0.01);
//    cg1.setTimeStep(0.004);
//    cg1.setTimeMax(0.2);
//    cg1.statFile.setSaveCount(5);
//    cg1.statFile.setName(problem.getName() + ".LucyO.TS.stat");
//    problem.cgHandler.copyAndAddObject(cg1);
//
//    TimeAveragedCG<CGCoordinates::O> cg2;
//    cg2.setTimeMin(0.6);
//    cg2.statFile.setSaveCount(50);
//    cg2.statFile.setName(problem.getName() + ".LucyO.TA.stat");
//    problem.cgHandler.copyAndAddObject(cg2);
//
    CG<CGCoordinates::Z> cgB;
    cgB.setNZ(200);
    cgB.setWidth(0.3);
    cgB.statFile.setSaveCount(20000);
    cgB.statFile.setName(problem.getName() + ".LucyZ.stat");
    problem.cgHandler.copyAndAddObject(cgB);
//
//    CG<CGFunctions::LinearZ> cgC;
//    cgC.setNZ(200);
//    cgC.setWidth(0.3);
//    cgC.statFile.setSaveCount(20000);
//    cgC.statFile.setName(problem.getName() + ".LinearZ.stat");
//    problem.cgHandler.copyAndAddObject(cgC);
//
//    CG<CGFunctions::HeavisideZ> cgD;
//    cgD.setNZ(200);
//    cgD.setWidth(0.3);
//    cgD.statFile.setSaveCount(20000);
//    cgD.statFile.setName(problem.getName() + ".HeavisideZ.stat");
//    problem.cgHandler.copyAndAddObject(cgD);
//
//    CG<CGFunctions::GaussZ> cgE;
//    cgE.setNZ(200);
//    cgE.setWidth(0.1);
//    cgE.statFile.setSaveCount(20000);
//    cgE.statFile.setName(problem.getName() + ".GaussZ.stat");
//    problem.cgHandler.copyAndAddObject(cgE);
//
//    CG<CGCoordinates::XZ> cgF;
//    cgF.setNX(30);
//    cgF.setNZ(30);
//    cgF.setWidth(0.45);
//    cgF.statFile.setSaveCount(20000);
//    cgF.statFile.setName(problem.getName() + ".LucyXZ.stat");
//    problem.cgHandler.copyAndAddObject(cgF);
//
//    CG<CGFunctions::LinearXZ> cgG;
//    cgG.setNX(30);
//    cgG.setNZ(30);
//    cgG.setWidth(0.45);
//    cgG.statFile.setSaveCount(20000);
//    cgG.statFile.setName(problem.getName() + ".LinearXZ.stat");
//    problem.cgHandler.copyAndAddObject(cgG);
//
//    CG<CGFunctions::HeavisideXZ> cgH;
//    cgH.setNX(30);
//    cgH.setNZ(30);
//    cgH.setWidth(0.45);
//    cgH.statFile.setSaveCount(20000);
//    cgH.statFile.setName(problem.getName() + ".HeavisideXZ.stat");
//    problem.cgHandler.copyAndAddObject(cgH);
//
//    CG<CGFunctions::GaussXZ> cgI;
//    cgI.setNX(30);
//    cgI.setNZ(30);
//    cgI.setWidth(0.15);
//    cgI.statFile.setSaveCount(20000);
//    cgI.statFile.setName(problem.getName() + ".GaussXZ.stat");
//    problem.cgHandler.copyAndAddObject(cgI);
//
//    CG<CGCoordinates::XYZ> cgJ;
//    cgJ.setNX(10);
//    cgJ.setNY(10);
//    cgJ.setNZ(10);
//    cgJ.setWidth(0.45);
//    cgJ.statFile.setSaveCount(20000);
//    cgJ.statFile.setName(problem.getName() + ".LucyXYZ.stat");
//    problem.cgHandler.copyAndAddObject(cgJ);
//
//    CG<CGFunctions::LinearXYZ> cgK;
//    cgK.setNX(10);
//    cgK.setNY(10);
//    cgK.setNZ(10);
//    cgK.setWidth(0.45);
//    cgK.statFile.setSaveCount(20000);
//    cgK.statFile.setName(problem.getName() + ".LinearXYZ.stat");
//    problem.cgHandler.copyAndAddObject(cgK);
//
//    CG<CGFunctions::HeavisideXYZ> cgL;
//    cgL.setNX(10);
//    cgL.setNY(10);
//    cgL.setNZ(10);
//    cgL.setWidth(0.45);
//    cgL.statFile.setSaveCount(20000);
//    cgL.statFile.setName(problem.getName() + ".HeavisideXYZ.stat");
//    problem.cgHandler.copyAndAddObject(cgL);
//
//    CG<CGFunctions::GaussXYZ> cgM;
//    cgM.setNX(10);
//    cgM.setNY(10);
//    cgM.setNZ(10);
//    cgM.setWidth(0.15);
//    cgM.statFile.setSaveCount(20000);
//    cgM.statFile.setName(problem.getName() + ".GaussXYZ.stat");
//    problem.cgHandler.copyAndAddObject(cgM);
//
//    CG<CGFunctions::GaussY> cg7;
//    cg7.setNY(200);
//    cg7.setWidth(0.15);
//    cg7.statFile.setSaveCount(20000);
//    cg7.statFile.setName(problem.getName() + ".GaussY.stat");
//    problem.cgHandler.copyAndAddObject(cg7);
//
//    CG<CGFunctions::GaussX> cg8;
//    cg8.setNX(200);
//    cg8.setWidth(0.15);
//    cg8.statFile.setSaveCount(20000);
//    cg8.statFile.setName(problem.getName() + ".GaussX.stat");
//    problem.cgHandler.copyAndAddObject(cg8);
//
//    CG<CGFunctions::GaussYZ> cg9;
//    cg9.setNY(30);
//    cg9.setNZ(30);
//    cg9.setWidth(0.15);
//    cg9.statFile.setSaveCount(20000);
//    cg9.statFile.setName(problem.getName() + ".GaussYZ.stat");
//    problem.cgHandler.copyAndAddObject(cg9);
//
//    CG<CGFunctions::GaussXY> cgA;
//    cgA.setNX(30);
//    cgA.setNY(30);
//    cgA.setWidth(0.15);
//    cgA.statFile.setSaveCount(20000);
//    cgA.statFile.setName(problem.getName() + ".GaussXY.stat");
//    problem.cgHandler.copyAndAddObject(cgA);

}

void testCGHandler(DPMBase& problem)
{
    Mdouble totalMass;
    Vec3D momentumBalance;
    
    std::cout << "Analysing #0: spatially averaged statistics rho(t) at t=1" << std::endl;
    // cg0 stores the last available time step for analysis
    CG<CGCoordinates::O>& cg0 = 
        *dynamic_cast<CG<CGCoordinates::O>*>(problem.cgHandler.getObject(0));
    cg0.writeAll(std::cout);
    
    //check total mass
    totalMass = cg0.getPoint(0).getDensity() - 0.8;
    momentumBalance = cg0.getPoint(0).getInteractionForceDensity() 
        + cg0.getPoint(0).getDensity() * problem.getGravity();
    std::cout << "Error in totalMass=" << totalMass << std::endl;
    std::cout << "Error in momentumBalance=" << Vec3D::getLength(momentumBalance) << std::endl;
    
    //case B: z-statistics
    CG<CGCoordinates::Z>& cgB = 
        *dynamic_cast<CG<CGCoordinates::Z>*>(problem.cgHandler.getObject(1));
    cgB.writeAll(std::cout); std::cout << std::endl;
    totalMass = 0.0;
    momentumBalance.setZero();
    for (auto& p: cgB.getPoints())
    {
        std::cout << "z=" << p.coordinates.getZ() << std::endl;
        totalMass += p.getDensity();
        momentumBalance += p.getInteractionForceDensity() 
        + p.getDensity() * problem.getGravity();
    }
    totalMass = totalMass/cgB.getPoints().size() - 0.8;
    momentumBalance /= cgB.getPoints().size();
    std::cout << "Error in totalMass=" << totalMass << std::endl;
    std::cout << "Error in momentumBalance=" << Vec3D::getLength(momentumBalance) << std::endl;
 
    
}

int main(int argc, char *argv[])
{
    //declare the DPM problem and set the name
    CGStaticBalanceSelfTest problem;
    problem.setName("CGStaticBalanceSelfTest");

    //problem is steady after t=0.6, see gnuplot:
    //set logscale y; p 'CGStaticBalanceSelfTest.ene' u 1:3
    setCGHandler(problem);
    
    //run the simulation
    problem.solve(argc, argv);
    
    //problem is steady after t=0.6, see gnuplot:
    //set logscale y; p 'CGStaticBalanceSelfTest.ene' u 1:3
    testCGHandler(problem);

}

