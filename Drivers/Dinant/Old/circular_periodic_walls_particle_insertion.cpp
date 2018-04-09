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

#include "scr/Mercury3D.h"
#include <sstream>
#include <iostream>
#include <cstdlib>

class circular_periodic_walls_particle_insertion: public Mercury3D {
  
	void setupInitialConditions()
	{
		particleHandler.clear();
		get_BoundaryHandler().clear();
		wallHandler.clear();
		
		CircularPeriodicBoundary b0;
		get_BoundaryHandler().copyAndAddBoundary(b0);
		
		CylindricalWall c0;
		c0.set(8);
		wallHandler.copyAndAddObject(c0);
		
		InfiniteWall w0;
		w0.set(Vec3D( 0, 0, -1), -getZMin());
		wallHandler.copyAndAddObject(w0);
		
		
		P0.setRadius(BaseParticleRadius);
		P0.computeMass(Species);
		P0.setPosition(Vec3D(0.0,0.0,getZMax()));
		P0.setVelocity(Vec3D(0.0,0.0,0.0));
		particleHandler.copyAndAddObject(P0);
	}
	
	protected:
	///Couts time
	virtual void printTime() const {
		cout << "t=" << setprecision(3) <<fixed<< left << setw(6) << getTime() 
			<< ", tmax=" << setprecision(3) << left << setw(6) << getTimeMax() << endl;
		cout.flush();
	}
	
	void actionsBeforeTimeStep(){
		
		int failed = 0;
		double createSize=0.4;
		while (failed<=10)
		{
			P0.setPosition(Vec3D(random.getRandomNumber(-createSize,createSize),random.getRandomNumber(-createSize,createSize),random.getRandomNumber(getZMax()-createSize,getZMax()+createSize)));
			if (InsertIfPossible(P0))
			{
				failed = 0; 
			}
			else failed++;
		}
	}	
	
	bool InsertIfPossible(BaseParticle &P)
	{
		particleHandler.copyAndAddObject(P);
		if(TestObjAgainstGrid(grid,particleHandler.back()))
		{
			//cout<<"Index of newest created Particle="<<particleHandler.back()->getIndex()<<endl;
			return true;
		}
		else
		{
			///todo{Maybe also check if the last particle in a Particular level is removed}
			Cell cell(particleHandler.back()->get_HGRID_x(),particleHandler.back()->get_HGRID_y(), particleHandler.back()->get_HGRID_z(), particleHandler.back()->get_HGRID_Level());
			int bucket = grid->ComputeHashBucketIndex(cell);
			grid->firstBaseParticleInBucket_[bucket] = particleHandler.back()->get_HGRID_NextObject();
			particleHandler.removeLastParticle();
			return false;
		}
	}
	Particle P0;
	
	public:
		double ParticleRadius;		
	
};

int main(int /*argc*/, char **/*argv[]*/)
{
	circular_periodic_walls_particle_insertion problem;
	problem.setName("circular_periodic_walls_particle_insertion");
	problem.setDensity(1000);
	problem.setXMin(-1);
	problem.setYMin(-2);
	problem.setXMax(8);
	problem.setYMax(8);
	problem.setZMin(0);
	problem.setZMax(4);
	problem.ParticleRadius=0.2;
	problem.species->setCollisionTimeAndRestitutionCoefficient(0.01,0.4, std::pow(problem.ParticleRadius,3)*constants::pi*4.0/3.0*1000);
	problem.setGravity(Vec3D(0.0,0.0,-10.0));
	problem.setTimeMax(2.5);
	problem.setTimeStep();
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimestep(1000,problem.getTimeMax(),problem.getTimeStep()));
	problem.write(std::cout,false);
	problem.solve();	
}
