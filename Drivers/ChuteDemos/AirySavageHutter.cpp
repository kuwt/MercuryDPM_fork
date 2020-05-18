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

#include <iostream>
#include "ChuteWithHopper.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"

///This code does the MD of a normal shock into a wall.
class AirySavageHutter : public ChuteWithHopper
{
public:
    
    void setupInitialConditions() override
    {
        logger(INFO, "Entering the solve now what happens");
        logger(INFO, "Problem name %", getName());
        
        //Setup the base i.e. the chute particles
        ChuteWithHopper::setupInitialConditions();
    }
    
    void actionsOnRestart() override
    {        
        logger.assert_always(getTimeMax() > getTime(), "Problem has been run and is complete: About to quit");
        logger(INFO, "Problem AirySavageHutter is not complete, will restart shortly. Current status:");        
        write(std::cout, false);
    }
    
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    AirySavageHutter problem;
    
    // Problem parameters
    problem.setName("AirySavageHutter");
    problem.setTimeStep(1e-4);
    problem.setTimeMax(500.0);
    problem.setHGridMaxLevels(2);
    
    // Particle properties
    problem.setInflowParticleRadius(0.5);
    problem.setFixedParticleRadius(0.25);
    problem.setRoughBottomType(MULTILAYER);
    LinearViscoelasticFrictionSpecies species;
    species.setDensity(6 / constants::pi);
    species.setStiffness(2e5);
    species.setSlidingStiffness(2.0 / 7.0 * species.getStiffness());
    species.setDissipation(25.0);
    species.setSlidingDissipation(0);
    species.setSlidingFrictionCoefficient(0.5);
    problem.speciesHandler.copyAndAddObject(species);
    
    // Chute properties
    problem.setChuteAngleAndMagnitudeOfGravity(27.0, 1.0);
    problem.setChuteLength(250);
    problem.setChuteWidth(10);
    problem.setMaxFailed(6);
    problem.makeChutePeriodic();
    
    //Hopper properties
    Mdouble exitHeight = 10.0;
    Mdouble exitLength = 10.0;
    Mdouble hopperAngle = 45.0;
    Mdouble hopperLength = 2.0 * exitLength;
    Mdouble hopperHeight = hopperLength;
    problem.setHopper(exitLength, exitHeight, hopperAngle, hopperLength, hopperHeight);
    problem.setHopperFillingPercentage(50.0);
    
    // Xballs tunning
    problem.setXBallsScale(0.012);
    problem.setXBallsAdditionalArguments("-sort -v0 -oh 500 -cmode 4");
    
    // You really do not need data more often that once a second do not change this number away from this without 
    // really good reason.
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(500, problem.getTimeMax(),
                                                                                     problem.getTimeStep()));
    logger(INFO, "time step = %", problem.getTimeStep());
    
    problem.solve();
}
