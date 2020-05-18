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

//! [RB:headers]
#include <sstream>
#include <Species/LinearViscoelasticSpecies.h>
#include "Chute.h"
//! [RB:headers]

/*!
 * This is a demo on how to implement a chute flow with a rough bottom. Note
 * that this is the same code as roughBottomSelfTest, but with a bigger
 * chute, a larger gravity and a longer simulation time.
 */

//! [RB:main]
int main()
{
    //Print description
    logger(INFO, "\nDemo of the chute flow with a rough bottom");
    
    //Construct the problem and assign a name
    Chute roughBottomSelfTest;
    roughBottomSelfTest.setName("roughBottomMultiLayer");


    //Set time stepping parameters
    roughBottomSelfTest.setTimeStep(1e-4);
    roughBottomSelfTest.setTimeMax(0.1);

    //Set output parameter: write to the output files every 50 time steps
    roughBottomSelfTest.setSaveCount(50);

    //Set the particle radii and the type of the rough bottom
    roughBottomSelfTest.setInflowParticleRadius(0.5);
    roughBottomSelfTest.setFixedParticleRadius(0.5);
    roughBottomSelfTest.setRoughBottomType(MULTILAYER);

    //Set the particle properties
    LinearViscoelasticSpecies* species = roughBottomSelfTest.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    species->setDensity(6.0/constants::pi);
    species->setStiffness(2e5);
    species->setDissipation(25.0);

    //Set the chute properties
    roughBottomSelfTest.setChuteAngleAndMagnitudeOfGravity(24.0, 1.0);
    roughBottomSelfTest.setChuteLength(3.0);
    roughBottomSelfTest.setChuteWidth(3.0);

    //Solve the system
    roughBottomSelfTest.solve();
}
//! [RB:main]
