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

#include "scr/Mercury3D.h"
#include "scr/StatisticsVector.h"
#include <sstream>
#include <iostream>
#include <cstdlib>

template <StatType T> class statistics_while_running : public StatisticsVector<T>, public Mercury3D
{
	public:
	statistics_while_running<T>() : StatisticsVector<T>() , Mercury3D() {}
	
	void setupInitialConditions() 
	{
		setXMax(1);
		setYMax(1.2);
		setZMax(1);


		set_N(3);
		particleHandler.getObject(0)->->setPosition(Vec3D(0.4, 0.4, 0.0));
		particleHandler.getObject(0)->->setVelocity(Vec3D(0.0, 1.0, 0.0));
		particleHandler.getObject(0)->->getRadius()=0.1;
		
		particleHandler.getObject(1)->setPosition(Vec3D(0.4, 0.8, 0.0));
		particleHandler.getObject(1)->setVelocity(Vec3D(0.0,-1.0, 0.0));
		particleHandler.getObject(1)->getRadius()=0.1;
		
		particleHandler.getObject(2)->Position=Vec3D(0.2, 0.7, 0.0);
		particleHandler.getObject(2)->Velocity=Vec3D(-1.0,-1.0, 0.0);
		particleHandler.getObject(2)->Radius=0.1;
		
				
		set_NWall(2);
		wallHandler.getObject(0)->set(Vec3D(0.0, 1.0,0.0), getYMax());
		Walls[1].set(Vec3D(0.0,-1.0,0.0),-getYMin());

		set_NWallPeriodic(1);
		WallsPeriodic[0].set(Vec3D(1,0,0), getXMin(), getXMax());
	}
};

int main(int argc, char* argv[])
{
	
    statistics_while_running<Y> problem;
    
    problem.setDoPeriodicWalls(false);
    //~ problem.doTimeAverage(false);
    problem.setN(20);
    problem.setTimeMaxStat(0.01);
    problem.setCGTimeMin(0.0);
	problem.setDensity(20);
	problem.setCGWidth(0.1);
	problem.setStiffnessAndRestitutionCoefficient(10000,0.8,problem.get_mass_from_Radius(0.0037));		
	problem.speciesHandler.getObject(0)->setSlidingFrictionCoefficient(0.5);
	problem.speciesHandler.getObject(0)->setSlidingStiffness(problem.speciesHandler.getObject(0)->getStiffness());
	problem.dataFile.setFileType(FileType::ONE_FILE);
	problem.fStatFile.setFileType(FileType::ONE_FILE);
	
	problem.setTimeMax(0.01);
	problem.setGravity(Vec3D(0,0,0));
	
	problem.setTimeStep();
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(3,problem.getTimeMax(),problem.getTimeStep()));
	problem.setName("statistics_while_running");		
	problem.write(std::cout,false);
	
	problem.readArguments(argc, argv);
	problem.solve();
}

