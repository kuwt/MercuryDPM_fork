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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ChuteWithHopper.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"


/**
 * A three-dimensional hopper inflow.
 **/
int main(int argc, char* argv[])
{
    //Print description
    logger(INFO, "Description: A three-dimensional hopper inflow.\n\n");
    
    // Problem parameters
    ChuteWithHopper problem;
    problem.setName("hopper_3d_demo");
    problem.setTimeMax(2.0); //Should be 10 for full length problem, but here I keep it low for a test case
    
    // Particle properties
    LinearViscoelasticSlidingFrictionSpecies species;
    species.setDensity(2400.0);
    species.setCollisionTimeAndRestitutionCoefficient(4e-3, 0.6, 2400);
    species.setSlidingDissipation(species.getDissipation());
    species.setSlidingFrictionCoefficient(0.8);
    problem.speciesHandler.copyAndAddObject(species);
    
    
    // Chute properties
    problem.setChuteAngle(30.0);
    problem.setChuteLength(50.0);
    problem.setChuteWidth(100);
    problem.setInflowParticleRadius(0.5, 0.5);
    problem.setFixedParticleRadius(0.0);
    problem.setRoughBottomType(MONOLAYER_ORDERED);
    
    problem.setMaxFailed(6);
    problem.makeChutePeriodic();
    problem.setTimeStep(4e-3 / 50);
    
    // Hopper properties
    problem.setHopperDimension(2);
    problem.setHopperLift(20);
    Mdouble ExitHeight = 4.0, ExitLength = 1.0 * ExitHeight, hopperAngle_ = 45.0, hopperLength_ = 3.0 * ExitLength;
    Mdouble hopperLowestPoint_ = ExitHeight - ExitLength * tan(problem.getChuteAngle());
    Mdouble hopperHeight_ =
            hopperLowestPoint_ + 1.1 * 0.5 * (hopperLength_ + ExitLength) / tan(hopperAngle_ * constants::pi / 180.0);
    Mdouble HopperCornerHeight =
            hopperHeight_ - 0.5 * (hopperLength_ - ExitLength) / tan(hopperAngle_ * constants::pi / 180.0);
    if (HopperCornerHeight <= 0.0)
    {
        hopperHeight_ += -HopperCornerHeight + problem.getInflowParticleRadius();
        HopperCornerHeight = problem.getInflowParticleRadius();
    }
    problem.setHopper(ExitLength, ExitHeight, hopperAngle_, hopperLength_, hopperHeight_);
    
    // Output properties
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(100, problem.getTimeMax(),
                                                                                     problem.getTimeStep()));
    problem.restartFile.setSaveCount(problem.getTimeMax() / problem.getTimeStep());
    problem.setXBallsAdditionalArguments("-sort -v0 -solidf -drotphi 0.05 -v0 -oh -200 -p 20 -noborder 3");
    problem.setParticlesWriteVTK(true);
    
    //solve
    problem.solve(argc, argv);
}
