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

/// In this file basic external statistics are checked by three vertically 
/// aligned particles, where the lowest particle is fixed and the other two 
/// are places such that the system is in equilibirum. XYZ and Z statistics 
/// are calculated for Gaussian, HeavisideSphere, Heaviside, Linear and Lucy
class myproblem : public DPMBase {
	
	void setupInitialConditions()
	{
		setSystemDimensions(3);
		setParticleDimensions(3);
		
		setXMax(2);
		setYMax(2);
		setZMax(2);
		setZMin(-1);
				
		double angle = 10.*constants::pi/180.;
		setGravity(Vec3D(0,sin(angle),-cos(angle)));
		setDensity(6/constants::pi);
		setCollisionTimeAndRestitutionCoefficient(1,sqrt(0.5),1);
		setSlidingFrictionCoefficient(0);
		setTimeStep(1e-2);
		setTimeMax(1);
		setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(3,getTimeMax(),getTimeStep()));
	
		//Create the particles	
		SphericalParticle p0;
	
		//overlap at which the system is in balance
		double delta = 1./getStiffness();
		
		//Particle 1 
		p0.setVelocity(Vec3D(0.0,0.0,0.0));
		p0.setRadius(0.5);	
		p0.setPosition(Vec3D(1,1,-0.5)-2.0*(p0.getRadius()-delta)*getGravity());
		particleHandler.copyAndAddObject(p0);
		
		//Particle 2
		p0.setPosition(Vec3D(1,1,-0.5)-(4.0*p0.getRadius()-3.0*delta)*getGravity());
		particleHandler.copyAndAddObject(p0);
		
		
		//Particle 3
		p0.setPosition(Vec3D(1,1,-0.5));
		p0.fixParticle();
		particleHandler.copyAndAddObject(p0);
	}
	
};

///check different CG
void testXYZ() {
	const StatType T=XYZ;
	int nz = 400, nx = 1, ny = 1;
	double w = 0.5; //w_over_rmax
	int verbosity = 0;
	
	std::cout << std::endl << "GAUSSIAN-XYZ---------------------" << std::endl;
	StatisticsVector<T> statsXYZ("check_statistics2");
	statsXYZ.getStatFile().setName("check_statistics2_XYZ_Gaussian.stat");
	statsXYZ.setVerbosityLevel(verbosity);
	statsXYZ.setNX(nx);
	statsXYZ.setNY(ny);
	statsXYZ.setNZ(nz);
	statsXYZ.setCGWidth_over_rmax(w/3.);
	statsXYZ.statistics_from_fstat_and_data();
	
	std::cout << std::endl << "Heaviside-XYZ-------------------------" << std::endl;
	StatisticsVector<T> statsXYZ_H("check_statistics2");
	statsXYZ_H.getStatFile().setName("check_statistics2_Heaviside.stat");
	statsXYZ_H.setVerbosityLevel(verbosity);
	statsXYZ_H.setCGShape(HeavisideSphere);
	statsXYZ_H.setNX(nx);
	statsXYZ_H.setNY(ny);
	statsXYZ_H.setNZ(nz);
	statsXYZ_H.setCGWidth_over_rmax(w);
	statsXYZ_H.statistics_from_fstat_and_data();
	
	std::cout << std::endl << "Heaviside_Polynomial-XYZ--------------" << std::endl;
	StatisticsVector<T> statsXYZ_P("check_statistics2");
	statsXYZ_P.getStatFile().setName("check_statistics2_XYZ_HeavisidePolynomial.stat");
	statsXYZ_P.setVerbosityLevel(verbosity);
	statsXYZ_P.setCGShape("Heaviside");
	statsXYZ_P.setNX(nx);
	statsXYZ_P.setNY(ny);
	statsXYZ_P.setNZ(nz);
	statsXYZ_P.setCGWidth_over_rmax(w);
	statsXYZ_P.statistics_from_fstat_and_data();
	
	std::cout << std::endl << "Linear_Polynomial-XYZ--------------" << std::endl;
	StatisticsVector<T> statsXYZ_l("check_statistics2");
	statsXYZ_l.getStatFile().setName("check_statistics2_XYZ_Linear.stat");
	statsXYZ_l.setVerbosityLevel(verbosity);
	statsXYZ_l.setCGShape("Linear");
	statsXYZ_l.setNX(nx);
	statsXYZ_l.setNY(ny);
	statsXYZ_l.setNZ(nz);
	statsXYZ_l.setCGWidth_over_rmax(w);
	statsXYZ_l.statistics_from_fstat_and_data();
	
	std::cout << std::endl << "Lucy-XYZ-------------------------" << std::endl;
	StatisticsVector<T> statsXYZ_L("check_statistics2");
	statsXYZ_L.getStatFile().setName("check_statistics2_XYZ_Lucy.stat");
	statsXYZ_L.setVerbosityLevel(verbosity);
	statsXYZ_L.setCGShape("Lucy");
	statsXYZ_L.setNX(nx);
	statsXYZ_L.setNY(ny);
	statsXYZ_L.setNZ(nz);
	statsXYZ_L.setCGWidth_over_rmax(w);
	statsXYZ_L.statistics_from_fstat_and_data();
}

///check different CG
void testZ() {
	const StatType T=Z;
	int nz = 400, nx = 1, ny = 1;
	double w = 0.5; //w_over_rmax
	int verbosity = 0;
	
	std::cout << std::endl << "GAUSSIAN-Z-------------------------" << std::endl;
	StatisticsVector<T> statsZ("check_statistics2");
	statsZ.getStatFile().setName("check_statistics2_Z_Gaussian.stat");
	statsZ.setVerbosityLevel(verbosity);
	statsZ.setNX(nx);
	statsZ.setNY(ny);
	statsZ.setNZ(nz);
	statsZ.setCGWidth_over_rmax(w/3.);
	statsZ.statistics_from_fstat_and_data();
	
	std::cout << std::endl << "Heaviside-Z-------------------------" << std::endl;
	StatisticsVector<T> statsZ_H("check_statistics2");
	statsZ_H.getStatFile().setName("check_statistics2_Z_Heaviside.stat");
	statsZ_H.setVerbosityLevel(verbosity);
	statsZ_H.setCGShape(HeavisideSphere);
	statsZ_H.setNX(nx);
	statsZ_H.setNY(ny);
	statsZ_H.setNZ(nz);
	statsZ_H.setCGWidth_over_rmax(w);
	statsZ_H.statistics_from_fstat_and_data();
	
	std::cout << std::endl << "Heaviside_Polynomial---------------" << std::endl;
	StatisticsVector<T> statsZ_P("check_statistics2");
	statsZ_P.getStatFile().setName("check_statistics2_Z_HeavisidePolynomial.stat");
	statsZ_P.setVerbosityLevel(verbosity);
	statsZ_P.setCGShape("Heaviside");
	statsZ_P.setNX(nx);
	statsZ_P.setNY(ny);
	statsZ_P.setNZ(nz);
	statsZ_P.setCGWidth_over_rmax(w);
	statsZ_P.statistics_from_fstat_and_data();
	
	std::cout << std::endl << "Linear_Polynomial---------------" << std::endl;
	StatisticsVector<T> statsZ_l("check_statistics2");
	statsZ_l.getStatFile().setName("check_statistics2_Z_Linear.stat");
	statsZ_l.setVerbosityLevel(verbosity);
	statsZ_l.setCGShape("Linear");
	statsZ_l.setNX(nx);
	statsZ_l.setNY(ny);
	statsZ_l.setNZ(nz);
	statsZ_l.setCGWidth_over_rmax(w);
	statsZ_l.statistics_from_fstat_and_data();
	
	std::cout << std::endl << "Lucy-Z-------------------------" << std::endl;
	StatisticsVector<T> statsZ_L("check_statistics2");
	statsZ_L.getStatFile().setName("check_statistics2_Z_Lucy.stat");
	statsZ_L.setVerbosityLevel(verbosity);
	statsZ_L.setCGShape("Lucy");
	statsZ_L.setNX(nx);
	statsZ_L.setNY(ny);
	statsZ_L.setNZ(nz);
	statsZ_L.setCGWidth_over_rmax(w);
	statsZ_L.statistics_from_fstat_and_data();
}

int main(int argc UNUSED, char *argv[] UNUSED)
{
	myproblem problem;
	problem.set_name("check_statistics2");
	problem.solve();
	
	testXYZ();
	testZ();
	
	std::cout << std::endl << "Exact--------------------------" << std::endl;
	double delta = 1./problem.speciesHandler.getObject(0)->getStiffness();
	std::cout << std::endl << "Exact values:" << std::setprecision(7)
	<< std::endl << "Nu =" << "constants::pi/36 =" << constants::pi/36
	<< std::endl << "Density =" << "2/12 =" << 2./12
	<< std::endl << "NormalStress =" << "-g*[(1-delta)+(1-2*delta)]/16 =" << -problem.getGravity()*((1-delta)+(1-2*delta))/12
	<< std::endl;
	
	
}
