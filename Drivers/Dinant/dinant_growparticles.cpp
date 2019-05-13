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

#include "DPMBase.h"
#include "Math/ExtendedMath.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include <sstream>
#include <iostream>
#include <cstdlib>
#include "Species/LinearViscoelasticSpecies.h"

class dinant_growparticles : public DPMBase
{
  public:
    dinant_growparticles()
    {
        species  = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    }

	void computeExternalForces(BaseParticle* CI)
	{
		/// Now add on gravity
		CI->addForce(getGravity() * CI->getMass());
		///Finally walls
		//computeForcesDueToWalls(CI);
		///Background friction
        CI->addForce (-CI->getVelocity()*0.1* species->getDissipation());
	}

	void setupInitialConditions() 
	{
		std::cout<<"Run Particle Initial Conditions"<<std::endl;
		setXMin(0);
		setXMax(0.1);
		setYMin(0);
		setYMax(0.1);
		particleHandler.clear();
        SphericalParticle P0;
		for(int i=0;i<100;i++)
		{
			P0.setPosition(Vec3D(random.getRandomNumber(getXMin(),getXMax()),random.getRandomNumber(getYMin(),getYMax()),0));
			P0.setVelocity(Vec3D(0.0,0.0,0.0));
			endradius[i]=random.getRandomNumber(0.0037,2.0*0.0037);
			P0.setRadius(0.1*endradius[i]);
		}
		expfact=30;
					
		wallHandler.clear();
        boundaryHandler.clear();
        if(periodic)
        {
            PeriodicBoundary b0;
            b0.set(Vec3D(1,0,0), getXMin(), getXMax());
            boundaryHandler.copyAndAddObject(b0);
            b0.set(Vec3D(0,1,0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(b0);
        }
        else
        {
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
	}

	
	protected:
	
	void actionsBeforeTimeStep()
	{
		double t= getTime();
		double tm=getTimeMax();
		double fac;
		if(2*t<tm)
			fac=0.1+0.9*(1.0-exp(-expfact*t*2.0/tm))/(1.0-exp(-expfact));
		else
			fac=1;
            
		for(unsigned int i=0;i<particleHandler.getNumberOfObjects();i++)
		{
			particleHandler.getObject(i)->setRadius(endradius[i]*fac);
			particleHandler.getObject(i)->getSpecies()->computeMass(particleHandler.getObject(i));
			hGridInsertParticle(particleHandler.getObject(i));
		}
		setTimeStep(0.02*helpers::computeCollisionTimeFromKAndDispAndEffectiveMass(species->getStiffness(), species->getDissipation(), 0.5*particleHandler.getSmallestParticle()->getMass()));
	}
		
	double expfact;
	double endradius[100];
	double fac;
	bool periodic;
public:
    LinearViscoelasticSpecies* species;
};

int main(int /*argc*/, char **/*argv[]*/) {
  dinant_growparticles problem;
	
	problem.species->setDensity(20);
	Mdouble mass = problem.species->getMassFromRadius(0.0037);
	problem.species->setStiffnessAndRestitutionCoefficient(10000,0.8,mass);
	problem.setTimeMax(2);
	problem.setGravity(Vec3D(0,0,0));
	problem.setSaveCount(100);
	problem.setName("dinant_growparticles");		
	//problem.setHGridMaxLevels(10);
	problem.write(std::cout,false);
	problem.solve();
	
}

