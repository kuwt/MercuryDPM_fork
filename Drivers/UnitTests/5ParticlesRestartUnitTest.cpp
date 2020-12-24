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

#include "DPMBase.h"
#include "Math/Helpers.h"

/// In this file, 5 Particles (4 fixed) are loaded from files 
/// "5Particles.ini" and "5Particles.restart". The particles are aligned 
/// such that the single nonfixed particle rotates sinusoidally without 
/// moving. This is to test the behaviour of the tangential spring and 
/// the file loading routines.
int main(int argc UNUSED, char *argv[] UNUSED)
{
	
	helpers::writeToFile("5Particles.ini",
		"5 0 0 0 0 1.6 1 1.6 \n"
		"0   0 0    0 0 0  0.5  0 0 0  0 0 0  0\n"
		"0   0 1.4  0 0 0  0.5  0 0 0  0 0 0  0\n"
		"1.4 0 0    0 0 0  0.5  0 0 0  0 0 0  0\n"
		"1.4 0 1.4  0 0 0  0.5  0 0 0  0 0 0  0\n"
		"0.7 0 0.7  0 0 0  0.5  0 0 0  0 30 0  0\n"
		);

	helpers::writeToFile("5Particles.restart",
        "MercuryDPM 0.14 name 5ParticlesRestartUnitTest_restart revision 4946 repository https://svn.mercurydpm.org/SourceCode/Trunk\n"
        "dataFile    fileType ONE_FILE saveCount 5 counter 202 lastSavedTimeStep 1001\n"
        "fStatFile   fileType ONE_FILE saveCount 5 counter 202 lastSavedTimeStep 1001\n"
        "eneFile     fileType ONE_FILE saveCount 1 counter 1002 lastSavedTimeStep 1001\n"
        "restartFile fileType ONE_FILE saveCount 5 counter 202 lastSavedTimeStep 1001\n"
        "statFile    fileType ONE_FILE saveCount 5 counter 0\n"
        "interactionFile fileType ONE_FILE saveCount 5 counter 1 lastSavedTimeStep 4294967295\n"
        "xMin 0 xMax 1.6 yMin 0 yMax 1 zMin 0 zMax 1.6\n"
        "timeStep 1e-05 time 0.0100099999999998 ntimeSteps 1001 timeMax 0.01\n"
        "systemDimensions 3 particleDimensions 3 gravity 0 0 0 writeVTK 0 NO_FILE NO_FILE 0 0 0 0 random  0 0 0 0 0 0 0 xBallsArguments  -v0 -solidf\n"
        "Species 1\n"
        "LinearViscoelasticSlidingFrictionSpecies id 0 density 1.9098593 stiffness 200000 dissipation 0 slidingStiffness 57142.857 slidingDissipation 0 slidingFrictionCoefficient 0.5 slidingFrictionCoefficientStatic 0.5\n"
        "Walls 0\n"
        "Boundaries 0\n"
        "Particles 5\n"
        "BaseParticle id 0 indSpecies 0 position -0.0376743764985181 0 -0.0178234366605382 orientation 0.999395051545515 0 -0.0347783114359788 0 velocity -4.35530922224593 0 -1.99883567233961 angularVelocity 0 -8.259265876579 0 force 0 0 0 torque 0 0 0 radius 0.5 invMass 1.00000000895498\n"
        "BaseParticle id 1 indSpecies 0 position -0.0178234366605382 0 1.43767437649853 orientation 0.999395051545515 0 -0.0347783114359777 0 velocity -1.99883567233959 0 4.35530922224583 angularVelocity 0 -8.25926587657869 0 force 0 0 0 torque 0 0 0 radius 0.5 invMass 1.00000000895498\n"
        "BaseParticle id 2 indSpecies 0 position 1.41782343666051 0 -0.0376743764985151 orientation 0.999395051545515 0 -0.0347783114359767 0 velocity 1.99883567233947 0 -4.35530922224564 angularVelocity 0 -8.25926587657844 0 force 0 0 0 torque 0 0 0 radius 0.5 invMass 1.00000000895498\n"
        "BaseParticle id 3 indSpecies 0 position 1.43767437649853 0 1.41782343666051 orientation 0.999395051545515 0 -0.034778311435975 0 velocity 4.35530922224547 0 1.99883567233942 angularVelocity 0 -8.25926587657801 0 force 0 0 0 torque 0 0 0 radius 0.5 invMass 1.00000000895498\n"
        "BaseParticle id 4 indSpecies 0 position 0.7 0 0.7 orientation 0.999939404927204 0 0.0110084728199945 0 velocity 5.94120193922054e-13 0 -1.04430653949196e-14 angularVelocity 0 -3.03706350631416 0 force 0 0 0 torque 0 0 0 radius 0.5 invMass 1.00000000895498\n"
        "Interactions 0"
		);
	
 	DPMBase problem;
 	problem.setName("5Particles");
    problem.readRestartFile();
	problem.readDataFile("5Particles.ini");
	problem.setAppend(false);
	problem.solve();
	
	//Now check the rotational energy in the system
	//gnuplot: p '5Particles.ene' u 1:4
    return 0;
}
