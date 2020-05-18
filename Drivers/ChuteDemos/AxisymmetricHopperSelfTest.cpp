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

#include "AxisymmetricHopper.h"

/** A three-dimensional, axisymmetric hopper.
 **/
int main(int argc, char *argv[])
{
    //Print description
    cout << endl << "Description: A quasi-2D inclined plane with inflow "
            "conditions on the left boundary, and deletion of particles when "
            "they exit the domain." << endl;
 
    // Problem parameters
    AxisymmetricHopper axisymmetricHopper;
    axisymmetricHopper.setName("AxisymmetricHopperSelfTest");
    axisymmetricHopper.setTimeMax(std::sqrt(6.0*axisymmetricHopper.getZMax())); //let the particles fall to at most 16 times the domain size

    // Particle properties
    axisymmetricHopper.setFixedParticleRadius(0);
    axisymmetricHopper.setInflowParticleRadius(0.5);

    // Chute properties
    axisymmetricHopper.setChuteAngleAndMagnitudeOfGravity(0.0,1.0);

    //Output parameters
    axisymmetricHopper.setLastSavedTimeStep(1);
    axisymmetricHopper.fStatFile.setFileType(FileType::NO_FILE);
    axisymmetricHopper.setXBallsAdditionalArguments("-sort -v0 -solidf");

    //solve
    axisymmetricHopper.solve(argc,argv);
}

