//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include "Chute.h"
using namespace std;

/** A quasi-2D inclined plane with inflow conditions on the left boundary, 
 * and deletion of particles when they exit the domain.
 **/
int main(int argc, char *argv[])
{
    //Print description
    logger(INFO, "\nDescription: A quasi-2D inclined plane with inflow conditions on the left boundary, and deletion of"
                 " particles when they exit the domain.");
    
    // Problem parameters
    Chute problem;
    problem.setName("chute_demo");
    problem.setTimeMax(0.5);
    
    // Particle properties
    problem.setFixedParticleRadius(0.001);
    problem.setInflowParticleRadius(0.001);
    
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    double mass = 0.5*species->getMassFromRadius(0.5*(problem.getMinInflowParticleRadius() + problem.getMaxInflowParticleRadius()));
    species->setCollisionTimeAndRestitutionCoefficient(2.5e-3,0.8,mass);
    species->setSlidingDissipation(0.0);
    problem.setTimeStep(0.02 * 2.5e-3);

    
    // Chute properties
    problem.setChuteAngle(30.0);
    problem.setXMax(0.1);
    problem.setYMax(2.0*problem.getInflowParticleRadius());
    
    // Inflow properties
    problem.setInflowHeight(0.1);
    problem.setInflowVelocity(0.1);
    problem.setInflowVelocityVariance(0.02);
    
    // Output properties
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(1, problem.getTimeMax(), problem.getTimeStep()));
    problem.dataFile.setSaveCount(100*problem.restartFile.getSaveCount());
    problem.eneFile.setSaveCount(100*problem.restartFile.getSaveCount());
//    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(1,problem.getTimeMax(),problem.getTimeStep())); //minimize output to the last time step
//    problem.set_number_of_saves_data(100); //allow enough data output so the evolution can be viewed in xballs
//    problem.set_number_of_saves_ene(100);

    //solve
    problem.solve(argc, argv);
}
