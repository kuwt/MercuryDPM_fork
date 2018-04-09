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

///todo{This code is not working as is wanted}
#include<iostream>
#include <Species/LinearViscoelasticSpecies.h>

#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "MercuryTime.h"

/*!
 * Simple contacts, without any walls, boundaries, output, hGrid
 */
class Contact : public Mercury3D
{
	Mdouble polydispersity_;

public:
	
	Contact (ParticleSpecies* s, Mdouble polydispersity) 
	{
		polydispersity_=polydispersity;
		speciesHandler.copyAndAddObject(s);
		setTimeStep(5e-5);
	    setTimeMax(1.0);
	    if (polydispersity==1.0)
	    	setName("SpeedTestContactsMonodisperse");
	    else
	    	setName("SpeedTestContactsPolydisperse"+std::to_string(polydispersity));
	}

	void computeExternalForces(BaseParticle *CI) override
	{
		DPMBase::computeExternalForces(CI);
		CI->addForce(-CI->getPosition());
	}	

	void setupInitialConditions() override
	{
	    unsigned N = 4;	//number of particles
		Mdouble r = 0.002;
		Mdouble rMax = r*sqrt(polydispersity_);
		Mdouble rMin = r/sqrt(polydispersity_);
		Mdouble h = std::max(2.0*rMax,3.0*r); //mesh size

	    setXMax(0.5*h*N);
	    setYMax(getXMax());
	    setZMax(getXMax());
	    setXMin(-getXMax());
	    setYMin(-getYMax());
	    setZMin(-getZMax());
	    setFileType(FileType::NO_FILE);
	    //eneFile.setFileType(FileType::ONE_FILE);
    	//dataFile.setFileType(FileType::ONE_FILE);
    	setSaveCount(100);
	    setDimension(3);
	    
		BaseParticle p0;		
        p0.setSpecies(speciesHandler.getObject(0));
		p0.setRadius(rMax);
		for (unsigned i = 0; i < N*N*N; i++)
		{
			unsigned ix = i%N;
			unsigned iz = i/N/N;
			unsigned iy = (i-ix-N*N*iz) / N;
			p0.setPosition(h*Vec3D(ix,iy,iz)-0.5*(N-1)*h*Vec3D(1,1,1));
			p0.setVelocity(0.01*Vec3D(random(-1,1),random(-1,1),random(-1,1)));
			particleHandler.copyAndAddObject(p0);
			p0.setRadius(random(rMin,rMax));
		}
	}

	void actionsAfterSolve() override
	{
		std::cout << "max diameter " << 2.0*particleHandler.getLargestParticle()->getRadius() << std::endl;
		std::cout << "min diameter " << 2.0*particleHandler.getSmallestParticle()->getRadius() << std::endl;
		hGridInfo();
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	std::cout <<
	"Non dissipative particles are simulated in a potential force field f_i=-r_i"
	"Thus, we get a clump of particles with a constant collision rate."
	"No walls, boundaries, output is added, testing the speed of particle collisions." 
	<< std::endl;
	
    Time time;

 	LinearViscoelasticSpecies s;
	s.setDensity(2000);
    s.setStiffness(1e3);
	//s.setDissipation(0.005);

    time.tic();
 	Contact mono(&s,1.0);
	mono.solve();    
    std::cout << "Total time to run monodisperse simulation: " << time.toc() << "s (Expected: 1.1s)" << std::endl;
    //expected time was measured on Thomas' mac

    time.tic();
 	Contact poly(&s,2.0);
	poly.setTimeMax(0.4);
	poly.solve();    
    std::cout << "Total time to run polydisperse simulation: " << time.toc() << "s (Expected: 1.1s)" << std::endl;

    time.tic();
 	Contact highPoly(&s,5.0);
 	highPoly.setTimeMax(0.3);
	highPoly.solve();    
    std::cout << "Total time to run highly polydisperse simulation: " << time.toc() << "s (Expected: 1.1s)" << std::endl;
}
