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

/** this file checks statistics for a chain of
 a) a wall at (0,0,-1)*x=0
 b) a small particle at r=1/2, z=1/2
 c) a large particle at r=1, z=2
 d) a large particle at r=1, z=4
 e) a fixed small particle at r=1/2, z=5.5
 Radii are increased by 1% to obtain forces between the particles; 
 spring constants are chosen (k=100) such that the force between two small particles is 1;
 (thus f=2 between large particles and f=1.5 between large and small particles, f=0.5 between wall and small particle) 
 density is kept such that the small particles have mass 1.
 Statistics are obtained in XYZ, XZ, and Z (*.XYZ.stat, ...). 
 IFD on small and large particles is evaluated (*.XYZs.stat, *.XYZl.stat, ...)
**/
class myMD : public DPMBase {
	
	void setupInitialConditions()
	{
		double f = 1.01; //factor by which particle radius is increases
		
		setSystemDimensions(3);
		setParticleDimensions(3);
		setDensity(6./constants::pi/mathsFunc::cubic(f));
		setStiffness(1./(f-1));
		setDissipation(0);
		setSlidingFrictionCoefficient(0);
		//setTimeStep_by_mass(8);
		setTimeStep(1e-100);
		setTimeMax(getTimeStep());
		setSaveCount(1e5);
		
		setXMin(-1);
		setYMin(-1);
		setZMin(0);
		setXMax(1);
		setYMax(1);
		setZMax(5.5);
		
		set_NWall(1);
		Walls[0].set(Vec3D(0,0,-1),0);
		set_NWallPeriodic(0);
		
		SphericalParticle P;

		P.setPosition(Vec3D(0,0,2));
		P.setVelocity(Vec3D(0,0,0));
		P.setRadius(1*f);
		particleHandler.copyAndAddObject(P);

		P.setPosition(Vec3D(0,0,4));
		P.setVelocity(Vec3D(1,0,0));
		P.setRadius(1*f);
		particleHandler.copyAndAddObject(P);

		P.setPosition(Vec3D(0,0,5.5));
		P.setVelocity(Vec3D(0,0,0));
		P.setRadius(.5*f);
		particleHandler.copyAndAddObject(P);

		P.setPosition(Vec3D(0,0,0.5));
		P.setVelocity(Vec3D(1,0,0));
		P.setRadius(.5*f);
		P.fixParticle();
		particleHandler.copyAndAddObject(P);

	}
	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	
	///Start off by solving the default md
	myMD md;
	md.set_name("bidisperse");
	md.solve();
	md.write(std::cout,false);
	
	StatisticsVector<Z> stats("bidisperse");
	stats.setN(500);
	stats.setCGWidth_over_rmax(0.01);
	stats.setZMinStat(stats.getZMin()-3.0*stats.getCGWidth_over_rmax());
	stats.setZMaxStat(stats.getZMax()+3.0*stats.getCGWidth_over_rmax());
	stats.setDoTimeAverage(false);
	stats.verbose();
	//stats.setStressTypeForFixedParticles(3);
	
	stats.getStatFile().setName("bidisperse.Z.stat");
	stats.statistics_from_fstat_and_data();
	
	stats.getStatFile().setName("bidisperse.Zl.stat");
	stats.set_rmin(0.75);
	stats.statistics_from_fstat_and_data();
	
	stats.getStatFile().setName("bidisperse.Zs.stat");
	stats.set_rmin(0);
	stats.set_rmax(0.75);
	stats.statistics_from_fstat_and_data();
	
	/* Gnuplot
	t='check_statistics_bidisperse/bidisperse.Z.stat'           
	tl='check_statistics_bidisperse/bidisperse.Zl.stat' 
	ts='check_statistics_bidisperse/bidisperse.Zs.stat' 
	Density=4                                          
	ContactStressZZ=41                                           
	TractionZ=51     

	X=1; Y=2; Z=3;
	Nu=4; Density=5;
	MomentumX=6; MomentumY=7; MomentumZ=8;
	NormalStressXX=33; NormalStressXY=34; NormalStressXZ=35; NormalStressYX=36; NormalStressYY=37; NormalStressYZ=38; NormalStressZX=39; NormalStressZY=40; NormalStressZZ=41;
	TangentialStressXX=42; TangentialStressXY=43; TangentialStressXZ=44; TangentialStressYX=45; TangentialStressYY=46; TangentialStressYZ=47; TangentialStressZX=48; TangentialStressZY=49; TangentialStressZZ=50;
	NormalTractionX=51; NormalTractionY=52; NormalTractionZ=53;
	TangentialTractionX=54; TangentialTractionY=55; TangentialTractionZ=56;
	FabricXX=57; FabricXY=58; FabricXZ=59; FabricYY=60; FabricYZ=61; FabricZZ=62;
	#ContactStressXX=($33+$42); ContactStressXY=($34+$43); ContactStressXZ=($35+$44); ContactStressYX=($36+$45); ContactStressYY=($37+$46); ContactStressYZ=($38+$47); ContactStressZX=($39+$48); ContactStressZY=($40+$49); ContactStressZZ=($41+$50);
	
	v=NormalStressZZ;
	p t u 3:v w lp, ts u 3:v w lp, tl u 3:v w lp
	*/
}
