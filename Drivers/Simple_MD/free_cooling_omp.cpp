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

#include "omp.h"
#include "DPMBase.h"
#include "Time.h"


///this tests if the force loop can be parallelized by a simple omp command


class serialMD : public DPMBase{
public:

    void setupInitialConditions() {
 	setDensity(6./constants::pi);
	setCollisionTimeAndRestitutionCoefficient(1,1,1);
	setTimeStep(.02);
 	setTimeMax(.5);
	setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(10,getTimeMax(),getTimeStep()));
	
	int N=mathsFunc::square(30);
	int N1=static_cast<int>(sqrt(N))+1;
	setXMax(N1);
	setYMax(N1);
	
	for (int i=0;i<N;i++) 
		{
	
	    int ix=static_cast<int>(i%N1);
	    int iy=static_cast<int>(i/N1);
    
	    double x=(ix+1);
	    double y=(iy+1);
    
    	BaseParticle p0;
    	
    	p0.setPosition(Vec3D(x+random.getRandomNumber(-.01,.01),y+random.getRandomNumber(-.01,.01),0.0));
		p0.setVelocity(Vec3D(0.0,0.0,0.0));
		p0.setRadius(0.51);
		particleHandler.copyAndAddObject(p0);
    
	    
		}	
    }

protected:
    void computeExternalForces(BaseParticle *P){computeWalls(P);}
	
};

class parallelMD : public serialMD{
public:

    void computeAllForces()
    {
	///Reset all forces to zero
	for (std::vector<BaseParticle*>::iterator it = particleHandler.begin(); it!=particleHandler.end(); ++it)
	{
		(*it)->setForce(Vec3D(0.0,0.0,0.0));
		(*it)->setTorque(Vec3D(0.0,0.0,0.0));
	} //end reset forces loop

	///Now loop over all particles contacts computing force contributions
	
	#pragma omp parallel for num_threads(4)
	for (unsigned int i=0;i<particleHandler.getNumberOfObjects();i++)	 
	//for (std::vector<Particle*>::iterator it = particleHandler.begin(); it!=particleHandler.end(); ++it)
	{
		BaseParticle* it=particleHandler.getObject(i);	
		///Now loop over all other particles looking for contacts
		computeInternalForces(it);
		//end inner loop over contacts.
			
		computeExternalForces(it);

	}
	
	
	//end outer loop over contacts.
    }

};


int main(int argc UNUSED, char *argv[] UNUSED)
{
    Time t;
    t.tic();
    serialMD problem;
    problem.random.setRandomSeed(0);
    problem.set_name("free_cooling_serial");
    problem.solve();
    std::cout << "Time " << t.toc() << std::endl;

    t.tic();
    parallelMD problem2;
    problem2.random.setRandomSeed(0);
    problem2.set_name("free_cooling_omp");
    problem2.solve();
    std::cout << "Time " << t.toc()*((double) 4.0) << " (parallel)" << std::endl;
    
    std::cout << system("diff *.fstat|wc") << std::endl;
}
