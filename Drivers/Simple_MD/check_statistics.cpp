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

#include "StatisticsVector.h"

/// In this file basic external statistics are checked by two particles who are in contact
class myproblem : public DPMBase {
	
	void setupInitialConditions()
	{
		setSystemDimensions(3);
		setParticleDimensions(3);

		setXMax(2);
		setYMax(1);
		setZMax(1);

		SphericalParticle P0,P1;

		P0.setPosition(Vec3D(0.51,0.0,0.5));
		P1.setPosition(Vec3D(1.49,0.0,0.5));
		
		P0.setVelocity(Vec3D(1.0,0.0,0.0));
		P1.setVelocity(Vec3D(1.0,0.0,0.0));

		
		P0.setRadius(.5);
		P1.setRadius(.5);
		
		particleHandler.copyAndAddObject(P0);
		particleHandler.copyAndAddObject(P1);
	}
	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	myproblem problem;
	problem.set_name("check_statistics");
	problem.setTimeStep(1e-10);
	problem.setTimeMax(1e-10);
	problem.solve();
	
	StatisticsVector<X> stats("check_statistics");
	stats.setN(100);
	stats.setCGWidth_over_rmax(1.0);
	stats.setDoTimeAverage(false);
	stats.setVerbosityLevel(0);
	stats.statistics_from_fstat_and_data();
	
}
