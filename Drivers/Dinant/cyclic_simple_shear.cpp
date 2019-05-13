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

#include "Mercury2D.h"
#include "StatisticsVector.h"
#include "Walls/InfiniteWall.h"
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "Species/LinearViscoelasticSpecies.h"

template <StatType T> class cyclic_simple_shear : public StatisticsVector<T>, public Mercury2D
{
	public:
	cyclic_simple_shear<T>() : StatisticsVector<T>() , Mercury2D()
    {
        species  = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    }

	void computeExternalForces(BaseParticle* CI)
	{
		/// Now add on gravity
		CI->addForce(getGravity() * CI->getMass());
		///Finally walls
		//computeForcesDueToWalls(CI);
		///Background friction
        CI->addForce (-CI->getVelocity()*0.1* species->getDissipation());
        CI->addTorque(-CI->getAngularVelocity()*0.1* species->getSlidingDissipation());
	}

	void setupInitialConditions() 
	{
		if(readDataFile("Startpos.dat",8))
			std::cout <<"inlezen gelukt"<< std::endl;
		else
		{
			std::cout <<"inlezen mislukt"<< std::endl;
			exit(1);
		}
				
		wallHandler.clear();       
        InfiniteWall w0;
        w0.set(Vec3D( 0.0, 1.0, 0.0), getYMax());
        wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0.0,-1.0, 0.0),getMin());
		wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D( 1.0, 0.0, 0.0), getXMax());
		wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(-1.0, 0.0, 0.0),getMin());
        wallHandler.copyAndAddObject(w0);
	}
    
	int cycles;
	
	protected:
	
	void actionsBeforeTimeStep()
	{
        double fact=atan(0.3);
        double alpha=fact*(1.0-cos(cycles*(getTime())/getTimeMax()*2*constants::pi));
        Vec3D n2=Vec3D( cos(alpha),-sin(alpha),0.0);
		Vec3D n3=Vec3D(-cos(alpha), sin(alpha),0.0);
	   	Vec3D p2=Vec3D(getXMax(),0.5*(getYMax()+getYMin()),0.0);
		Vec3D p3=Vec3D(getXMin(),0.5*(getYMax()+getYMin()),0.0); 	
        dynamic_cast<InfiniteWall*>(wallHandler.getObject(2))->set(n2,Vec3D::dot(n2,p2));
		dynamic_cast<InfiniteWall*>(wallHandler.getObject(3))->set(n3,Vec3D::dot(n3,p3));
	}
public:
    LinearViscoelasticSlidingFrictionSpecies* species;
};

int main(int /*argc*/, char **/*argv[]*/)
{
    cyclic_simple_shear<Z> problem;

    problem.cycles=2;
    problem.autoNumber();
    problem.setDoPeriodicWalls(false);
    problem.setDoTimeAverage(false);
     
    problem.species->setDensity(5);
	Mdouble mass = problem.species->getMassFromRadius(0.0064);
	problem.species->setStiffnessAndRestitutionCoefficient(1e6,0.8,mass);
    //problem.species->setStiffnessAndRestitutionCoefficient(1e4,0.8,problem.getDensity()*constants::pi*pow(0.0037,2));
	problem.species->setSlidingFrictionCoefficient(0.5);
	problem.species->setSlidingStiffness(problem.species->getStiffness()*2/7);
	problem.dataFile.setFileType(FileType::ONE_FILE);
	problem.fStatFile.setFileType(FileType::ONE_FILE);
	problem.setGravity(Vec3D(0,0,0));

    problem.setTimeStep(1e-6);
    //problem.setTimeStepByParticle(0.02*constants::pi/sqrt(2*problem.getStiffness()/problem.particleHandler.get_LightestParticle()->getMass()));
	problem.setTimeMax(problem.cycles*1e6* problem.getTimeStep()); //1e6 seem quite good
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(problem.cycles*10+1, problem.getTimeMax(), problem.getTimeStep()));
    problem.eneFile.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(problem.cycles*400+1, problem.getTimeMax(), problem.getTimeStep()));
	//problem.set_number_of_saves(problem.cycles*10+1);
	//problem.set_number_of_saves_ene(problem.cycles*400+1);
	problem.setName("cyclic_simple_shear");		
	problem.write(std::cout,false);
	problem.solve();
}

