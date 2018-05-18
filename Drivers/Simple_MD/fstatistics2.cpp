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
		switch(0) {
		case 0:
			//case 0: two particles touching each other
			setSystemDimensions(3);
			setParticleDimensions(3);
			
			setXMax(2);
			setYMax(1);
			setZMax(1);
			
			set_NWall(0);
			set_NWallPeriodic(0);
			
			set_N(2);
			particleHandler.getObject(0)->Position=Vec3D(0.51,0.5,0.5);
			particleHandler.getObject(1)->Position=Vec3D(1.49,0.5,0.5);
			
			particleHandler.getObject(0)->setVelocity(Vec3D(1.0,0.0,0.0));
			particleHandler.getObject(1)->setVelocity(Vec3D(1.0,0.0,0.0));
			
			particleHandler.getObject(0)->getRadius()=.5;
			particleHandler.getObject(1)->getRadius()=.5;

			break;
		case 1:
			set_N(3);
			
			particleHandler.getObject(1)->Position=Vec3D(7.5e-3,2e-3,0e-3);
			particleHandler.getObject(0)->Position=Vec3D(7.5e-3,8e-3,0e-3);
			particleHandler.getObject(2)->Position=Vec3D(2.4e-3,5e-3,0e-3);
			
			particleHandler.getObject(0)->setVelocity(Vec3D(-0.1,0.0,0.0));
			particleHandler.getObject(1)->setVelocity(Vec3D( 0.1,0.0,0.0));
			particleHandler.getObject(2)->Velocity=Vec3D( 0.0,0.0,0.0);
			
			particleHandler.getObject(0)->getRadius()=2.5e-3;
			particleHandler.getObject(1)->getRadius()=2.5e-3;
			particleHandler.getObject(2)->Radius=2.5e-3;

			set_NWallPeriodic(1);
			WallsPeriodic[0].set(Vec3D(0,1,0), getYMin(), getYMax());

			set_NWall(1);
			Walls[0].set(Vec3D(-1.0, 0.0, 0.0), -getXMin());
			break;
		case 2:
		{
			int n = 4;
			set_N(mathsFunc::square(n));
			
			for (int i=0; i<n; i++)
			for (int j=0; j<n; j++) {
				Particles[i*n+j].Position=Vec3D((0.5+i)*getXMax()/n,(0.0+j)*getYMax()/n,0.0);
				Particles[i*n+j].Velocity=Vec3D(0.0,0.0,0.0);
				Particles[i*n+j].Radius=2.5e-3;
			}
			
			set_NWallPeriodic(2);
			WallsPeriodic[0].set(Vec3D(1,0,0), getXMin(), getXMax());
			WallsPeriodic[1].set(Vec3D(0,1,0), getYMin(), getYMax());
			
			set_NWall(0);
			break;
		}
		case 3:
		
			setSystemDimensions(3);
			setParticleDimensions(3);
			
			setXMax(4e-3);
			setYMax(2e-3);
			setZMax(8e-3);
			
			set_NWall(1);
			Walls[0].set(Vec3D(-1.0, 0.0, 0.0), -getXMin());

			set_NWallPeriodic(1);
			WallsPeriodic[0].set(Vec3D(0,0,1), getZMin(), getZMax());

			set_N(7);
			
			for (int i=0; i<get_N(); i++) {
				particleHandler.getObject(i)->setVelocity(Vec3D(0.0,0.0,0.0));
				particleHandler.getObject(i)->getRadius()=1e-3;
			}
			
			particleHandler.getObject(0)->Position=Vec3D( .9e-3,1.0000001e-3,1.5e-3);
			particleHandler.getObject(1)->Position=Vec3D(  1e-3,1e-3,3.5e-3);
			particleHandler.getObject(2)->Position=Vec3D(2.8e-3,1e-3,3.5e-3);
			Particles[3].Position=Vec3D(  1e-3,1e-3,5.5e-3);
			Particles[4].Position=Vec3D(2.8e-3,1e-3,5.5e-3);
			Particles[5].Position=Vec3D(  1e-3,1e-3,7.5e-3);
			Particles[6].Position=Vec3D(2.8e-3,1e-3,7.5e-3);
			
			Particles[4].fixParticle();

			break;
		case 4:
		
			setSystemDimensions(2);
			setParticleDimensions(2);
			
			setXMax(4e-3);
			setYMax(8e-3);
			setZMax(2e-3);
			
			set_NWall(1);
			Walls[0].set(Vec3D(-1.0, 0.0, 0.0), -getXMin());

			set_NWallPeriodic(1);
			WallsPeriodic[0].set(Vec3D(0,1,0), getYMin(), getYMax());

			set_N(7);
			
			for (int i=0; i<get_N(); i++) {
				particleHandler.getObject(i)->getVelocity() = Vec3D(0.0,i,0.0);
				particleHandler.getObject(i)->getRadius()=.999999e-3;
			}
			
			particleHandler.getObject(0)->Position=Vec3D( .9e-3,1.5e-3,0e-3);
			particleHandler.getObject(1)->Position=Vec3D(  1e-3,3.5e-3,0e-3);
			particleHandler.getObject(2)->Position=Vec3D(2.8e-3,3.5e-3,0e-3);
			Particles[3].Position=Vec3D(  1e-3,5.5e-3,0e-3);
			Particles[4].Position=Vec3D(2.8e-3,5.5e-3,0e-3);
			Particles[5].Position=Vec3D(  1e-3,7.5e-3,0e-3);
			Particles[6].Position=Vec3D(2.8e-3,7.5e-3,0e-3);
			
			Particles[4].fixParticle();

			break;
		}
	}
	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	
	///Start off by solving the default problem
	myproblem problem;
	problem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
	problem.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
	problem.getStatFile().setFileType(FileType::MULTIPLE_FILES_PADDED);
	problem.set_name("elastic_collision");
	problem.setTimeStep(1e-10);
	problem.setSaveCount(3);
	problem.setTimeMax(9e-10);
	problem.solve();
	problem.writeRestartFile();
	problem.write(std::cout,false);
	
	StatisticsVector<XYZ> stats("elastic_collision");
	stats.set_h(.1);
	stats.setCGWidth(.5);
	stats.doTimeAverage();
	stats.statistics_from_fstat_and_data();
	
}
