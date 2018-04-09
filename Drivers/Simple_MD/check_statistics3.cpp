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

#include "scr/StatisticsVector.h"

class myproblem : public DPMBase {
	
	void setupInitialConditions()
	{
		setSystemDimensions(3);
		setParticleDimensions(3);
		
		setXMax(2);
		setYMax(2);
		setZMax(3);
		
		set_NWall(0);
		//~ Walls[0].set(Vec3D(0,0,-1),0);
		
		set_NWallPeriodic(0);
		
		setGravity(Vec3D(0,0,0));
		setDensity(6/constants::pi);
		double tc = 1e-8;
		setCollisionTimeAndRestitutionCoefficient(tc,sqrt(0.5),1);
		setSlidingFrictionCoefficient(0);
		setTimeStep(tc/50.);
		setTimeMax(120*getTimeStep());
		setSaveCount(60);
		
		set_N(2);
		particleHandler.getObject(0)->getRadius()=.5;		
		particleHandler.getObject(0)->setVelocity(Vec3D(-0.1,0.1,-1.0));
		particleHandler.getObject(0)->Position=Vec3D(1,1,2);
		particleHandler.getObject(1)->getRadius()=.5;		
		particleHandler.getObject(1)->setVelocity(Vec3D(0.0,0.0,0.0));
		particleHandler.getObject(1)->Position=Vec3D(1,1,1);
	}
	
};

///check different CG
void testZ() {
	const StatType T=Z;
	int nz = 500, nx = 1, ny = 1;
	double w = 0.6; //w_over_rmax
	int verbosity = 1;

	
	std::cout << std::endl << "GAUSSIAN-Z-------------------------" << std::endl;
	StatisticsVector<T> statsZ("check_statistics3");
	statsZ.getStatFile().setName("check_statistics3.stat");
	statsZ.setVerbosityLevel(verbosity);
	statsZ.setNX(nx);
	statsZ.setNY(ny);
	statsZ.setNZ(nz);
	statsZ.setCGWidth_over_rmax(w);
	statsZ.setSuperExact(true);
	statsZ.statistics_from_fstat_and_data();

	
	std::cout << std::endl << "GAUSSIAN-Z-Gradient-------------" << std::endl;
	StatisticsVector<T> statsZG("check_statistics3");
	statsZG.getStatFile().setName("check_statistics3_Gradient.stat");
	statsZG.setVerbosityLevel(verbosity);
	statsZG.setNX(nx);
	statsZG.setNY(ny);
	statsZG.setNZ(nz);
	statsZG.setCGWidth_over_rmax(w);
	statsZG.setSuperExact(true);
	statsZG.setDoGradient(true);
	statsZG.statistics_from_fstat_and_data();
}

int main(int argc UNUSED, char *argv[] UNUSED)
{
	myproblem problem;
	problem.set_name("check_statistics3");
	problem.solve();
	problem.write(std::cout, true);
	
	testZ();
	//~ data=loadstatistics('check_statistics3_Gradient.stat')
	//~ plot(data.z,data.StrainZZ./data.Density.^2)
	//~ plot(data.z,(data.DisplacementMomentumZ_dz.*data.Density-data.DisplacementMomentumZ.*data.Density_dz)./data.Density.^2)
	//~ plot(data.z,data.DisplacementMomentumZ./data.Density)
	//~ plot(data.z,data.MomentumZ./data.Density)
	
}
