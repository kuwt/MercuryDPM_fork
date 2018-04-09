//Copyright (c) 2013-2017, The MercuryDPM Developers Team. All rights reserved.
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

#include "Booster.h"

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble bladeInclination = 20.0*constants::pi/180.0; // default angle of 5.71059314 deg such that tan(angle)=0.1 (it's easy to check the normals then)
    Mdouble bladeThickness = 7.2e-2;
    Mdouble bladeWidth = 40.0/180.0*constants::pi*64.25e-2; //set in degrees
    Mdouble bladeHeight = 0.3*64.25e-2;
    unsigned int bladeNumber = 2;

    BoosterWithSmallAngledBlade booster(bladeThickness, bladeHeight, bladeWidth, bladeInclination, bladeNumber);

    booster.autoNumber();
    std::vector<int> studyNum = booster.get2DParametersFromRunNumber(3,3);

    //booster.bladeInclination_ = (10.0 + 10.0*(studyNum[2]-1.))*constants::pi/180.0;
    booster.bladeHeight_ = (0.15 + 0.15*(studyNum[2]-1.))*64.25e-2;
    booster.bladeNumber_ = studyNum[1];
    //booster.bladeWidth_ = ((20.0 + 20.0*(studyNum[1]-1.))/180.0)*constants::pi*64.25e-2;

    booster.drumInclination = 5.0*constants::pi/180.0;
    booster.revolutionsPerSecond = 8.0/60.0;
    booster.particleDiameter = 7.2e-2/1.2599;// cuberoot(2)=1.2599
    booster.particleSpecies->setDensity(2000.0);
    booster.collisionTime = 0.01; //soft particles, but relative overlap due to gravity, o/d ~ (tc/tg)^2 = 0.0136, is still small
    booster.restitutionCoefficient = 0.001; // very dissipative particles
    booster.setXBallsAdditionalArguments(" -v0 -solidf  -w 1200 -s 0.76 -noborder 3");
    booster.setTimeMax(5.0*60.0/8.0); //run 5 revolutions
    booster.solve();
    return 0;
}
