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

//! [CH:headers]
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "ChuteWithHopper.h"
//! [CH:headers]

// A quasi-2D inclined plane with hopper inflow conditions, 
// and deletion of particles when they exit the domain.

//! [CH:main]
int main() 
{
   
    //Problem parameters
    ChuteWithHopper problem;
    problem.setName("HopperSelfTest");
    Mdouble tc = 4.0e-4;                // particle collision time
    problem.setTimeStep(0.02 * tc);  
    problem.setTimeMax(0.01);


    //Particle properties
    problem.setInflowParticleRadius(0.32e-3, 0.40e-3);
    problem.setFixedParticleRadius(0.36e-3);

    LinearViscoelasticSlidingFrictionSpecies species;
    species.setHandler(&problem.speciesHandler);
    species.setDensity(2400.0);
    Mdouble meanRadius = 0.5 * (problem.getMinInflowParticleRadius() + problem.getMaxInflowParticleRadius());
    Mdouble mass = species.getMassFromRadius(meanRadius);
    species.setCollisionTimeAndRestitutionCoefficient(tc, 0.95, mass);
    species.setSlidingFrictionCoefficient(0.5);
    problem.speciesHandler.copyAndAddObject(species);

    //Chute properties
    problem.setRoughBottomType(MONOLAYER_ORDERED); // try also MONOLAYER_DISORDERED, MULTILAYER or FLAT
    problem.setChuteAngleAndMagnitudeOfGravity(10.0,9.8); // in degrees (relative to horizontal)
    problem.setChuteLength(40.0e-3);
    problem.setChuteWidth(3.2e-3 / 2.0);

    //Hopper properties
    Mdouble ExitHeight = 10.0e-3;
    Mdouble ExitLength = ExitHeight;
    Mdouble hopperAngle_ = 60.0;
    Mdouble hopperLength_ = 4.0 * ExitLength;
    
//    Mdouble hopperLowestPoint_ = ExitHeight - ExitLength * std::tan(problem.getChuteAngle());
//    Mdouble hopperHeight_ = hopperLowestPoint_ + 1.1 * 0.5 * (hopperLength_ + ExitLength) / std::tan(hopperAngle_ * constants::pi / 180.0); // NB: this is not DEDUCED, but SET!!! (BvdH)
    
    Mdouble hopperHeight_ = 24.1139e-3;
    problem.setHopper(ExitLength, ExitHeight, hopperAngle_, hopperLength_, hopperHeight_);

    // Output properties
    problem.setSaveCount(51);

    //'time step ratio': size of minimum particle over maximum distance travelled 
    // by any particle per time step (i.e., should be >> 1 )
    std::cout << "time step ratio: " << problem.getTimeStepRatio() << std::endl;
    
    //solve
    problem.solve();
}
//! [CH:main]
