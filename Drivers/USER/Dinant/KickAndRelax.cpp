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

#include "Mercury3D.h"
#include "StatisticsVector.h"
#include "Walls/InfiniteWall.h"
#include <sstream>
#include <iostream>
#include <cstdlib>
#include "Species/LinearViscoelasticSpecies.h"

class KickAndRelax:  public Mercury3D
{
public:

    KickAndRelax()
    {
        species  = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    }

	void setupInitialConditions() 
	{		
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
            w0.set(Vec3D( 0.0, 1.0, 0.0),getMax());
            wallHandler.copyAndAddObject(w0);
            w0.set(Vec3D( 0.0,-1.0, 0.0),getMin());
            wallHandler.copyAndAddObject(w0);
            w0.set(Vec3D( 1.0, 0.0, 0.0),getMax());
            wallHandler.copyAndAddObject(w0);
            w0.set(Vec3D(-1.0, 0.0, 0.0),getMin());
            wallHandler.copyAndAddObject(w0);
        }
        
        double xMomentum=0.0,newXMomentum=0.0;;
        double yMomentum=0.0,newYMomentum=0.0;;
        
        for (std::vector<BaseParticle*>::const_iterator it=particleHandler.begin();it!=particleHandler.end();it++)
        {
            (*it)->setVelocity(Vec3D(random.getRandomNumber(-maxInitialVelocity_,maxInitialVelocity_),random.getRandomNumber(-maxInitialVelocity_,maxInitialVelocity_),0.0));
            xMomentum+=(*it)->getVelocity().X*(*it)->getMass();
            yMomentum+=(*it)->getVelocity().Y*(*it)->getMass();
        }
        std::cout<<"Voor correctie xMomentum="<<xMomentum<<std::endl;
        std::cout<<"Voor correctie yMomentum="<<yMomentum<<std::endl;
        for (std::vector<BaseParticle*>::const_iterator it=particleHandler.begin();it!=particleHandler.end();it++)
        {
            Vec3D Vel=(*it)->getVelocity();
            Vel.X-=xMomentum/particleHandler.getNumberOfObjects()/(*it)->getMass();
            Vel.Y-=yMomentum/particleHandler.getNumberOfObjects()/(*it)->getMass();
            (*it)->setVelocity(Vel);
            newXMomentum+=(*it)->getVelocity().X*(*it)->getMass();
            newYMomentum+=(*it)->getVelocity().Y*(*it)->getMass();
        }
        std::cout<<"Na correctie xMomentum="<<newXMomentum<<std::endl;
        std::cout<<"Na correctie yMomentum="<<newYMomentum<<std::endl;
	}
    
	
	double maxInitialVelocity_;
	bool periodic;
public:
    LinearViscoelasticSpecies* species;
};

int main(int argc, char* argv[])
{
    KickAndRelax problem;
    
	if (argc>1) {
		problem.maxInitialVelocity_=atof(argv[1]);
        std::cout << "Using Maximum initial velocity "<<problem.maxInitialVelocity_<<std::endl;
	} else {
		std::cout << "Maximum initial velocity" << std::endl;
		exit(-1);
	}   
    
    problem.species->setDensity(100);
    problem.setParticleDimensions(2);
    if(problem.readDataFile("Startpos.dat",14)) //7 to read from minEpot, 14 from restart data
        std::cout <<"inlezen gelukt"<< std::endl;
    else
    {
        std::cout <<"inlezen mislukt"<< std::endl;
        exit(1);
    }

    double collisionTime_=1e-3;
    double coefficientOfRestitution_=0.8;
    double minParticleRadius=problem.particleHandler.getSmallestParticle()->getRadius();
    std::cout<<"minParticleRadius="<<minParticleRadius<<std::endl;
    Mdouble smallestMass = problem.species->getMassFromRadius(minParticleRadius);
    problem.species->setCollisionTimeAndRestitutionCoefficient(collisionTime_,coefficientOfRestitution_,smallestMass);
    problem.setTimeStep(problem.species->getCollisionTime(smallestMass));
    problem.setGravity(Vec3D(0,0,0));
	problem.setTimeMax(5e4* problem.getTimeStep());
    
	problem.setSaveCount(50);
	problem.eneFile.setSaveCount(5);
    problem.fStatFile.setSaveCount(1e9);
	problem.setName("KickAndRelax");		
	problem.write(std::cout,false);
	problem.solve();
    return 0;
}

