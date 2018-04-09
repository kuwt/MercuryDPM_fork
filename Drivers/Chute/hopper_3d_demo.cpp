//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ChuteWithHopper.h"

using namespace std;

/** A three-dimensional hopper inflow.
 **/
int main(int argc, char *argv[])
{
    //Print description
    cout << endl << "Description: A three-dimensional hopper inflow." << endl;

    // Problem parameters
    ChuteWithHopper problem;
    problem.setName("hopper_3d_demo");
    problem.setTimeMax(2.0); //Should be 10 for full length problem, but here I keep it low for a test case
        
    // Particle properties
    problem.setDensity(2400.0);
    problem.setInflowParticleRadius(0.5,0.5);
    problem.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(4e-3, 0.6);
    problem.setTimeStep(4e-3/10);
    problem.speciesHandler.getObject(0)->setSlidingDissipation(problem.get_dissipation());
    problem.speciesHandler.getObject(0)->setSlidingFrictionCoefficient(0.8);
    problem.setFixedParticleRadius(0.0);
    problem.setRoughBottomType(MONOLAYER_ORDERED);
    problem.setAlignBase(false);
    
    // Chute properties
    problem.setChuteAngle(30.0);
    problem.setChuteLength(50.0);
    problem.setChuteWidth(100);
    problem.setMaxFailed(6);
    problem.makeChutePeriodic();
    
    // Hopper properties
    problem.setHopperDimension(2);
    problem.setHopperLift(20);
    Mdouble ExitHeight = 4.0, ExitLength = 1.0 * ExitHeight, hopperAngle_ = 45.0, hopperLength_ = 3.0 * ExitLength;
    Mdouble hopperLowestPoint_ = ExitHeight - ExitLength * tan(problem.getChuteAngle());
    Mdouble hopperHeight_=hopperLowestPoint_ + 1.1 * 0.5*(hopperLength_+ExitLength) / tan(hopperAngle_*constants::pi/180.0);
    Mdouble HopperCornerHeight = hopperHeight_ - 0.5*(hopperLength_-ExitLength) / tan(hopperAngle_*constants::pi/180.0);
    if (HopperCornerHeight<=0.0) { hopperHeight_ += -HopperCornerHeight + problem.getInflowParticle()->getRadius(); HopperCornerHeight = problem.getInflowParticle()->getRadius(); }    
    problem.set_Hopper(ExitLength,ExitHeight,hopperAngle_,hopperLength_,hopperHeight_);    
        
    // Output properties
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimestep(1,problem.getTimeMax(),problem.getTimeStep())); //minimize output to the last timestep
    problem.set_number_of_saves_data(100); //allow enough data output so the evolution can be viewed in xballs
    problem.set_number_of_saves_ene(100);
    problem.setXBallsAdditionalArguments("-sort -v0 -solidf -drotphi 0.05 -v0 -oh -200 -p 20 -noborder 3");

    cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << endl; // speed allowed before particles move through each other!
    cout << "dt=" << problem.getTimeStep() << endl;

    //solve
    problem.solve(argc,argv);
}
