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

#include<iostream>
#include "Species/LinearViscoelasticSpecies.h"
#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "MercuryTime.h"

/*!
 * Simple contacts, without any walls, boundaries, output, hGrid
 */
class Contact : public Mercury3D
{
	Mdouble polydispersity_;

public:

    explicit Contact (Mdouble polydispersity)
	{
		LinearViscoelasticSpecies s;
		s.setDensity(2000);
		s.setStiffness(1e3);
		speciesHandler.copyAndAddObject(s);

		polydispersity_=polydispersity;
		setTimeStep(5e-5);
	    setTimeMax(0.39);
    	setName("SpeedTest_P"+helpers::to_string(polydispersity,2));
		setFileType(FileType::NO_FILE);
	}

	void computeExternalForces(BaseParticle *CI) override
	{
		DPMBase::computeExternalForces(CI);
		CI->addForce(-CI->getPosition());
	}	

	void setupInitialConditions() override
	{
	    unsigned N = 10;	//number of particles
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
    	setSaveCount(100);
	    setDimension(3);
	    
		SphericalParticle p0;
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
		logger(INFO,"%<r<%, N=%, C=%",
			   particleHandler.getSmallestInteractionRadius(),
			   particleHandler.getLargestInteractionRadius(),
			   particleHandler.getNumberOfObjects(),
			   interactionHandler.getNumberOfObjects()
		);
		//hGridInfo();
	}

};

int main()
{
	std::cout <<
	"A gas of non-dissipative particles are simulated, colliding at a constant rate.\n"
	"An external potential force field, f_i=-r_i, keeps the particles in the center of the box.\n"
 	"No walls, boundaries, file output, thus testing the speed of particle collisions."
	<< std::endl;
	
    Time time;

    time.tic();
 	Contact mono(1.0);
	mono.solve();
    std::cout << "Total time to run monodisperse simulation: " << time.toc() << "s (Expected: 3s)" << std::endl;
    //expected time was measured on Thomas' pc 26-Apr-2018 (r2816, Release)

    time.tic();
 	Contact poly(2.0);
	poly.setTimeMax(0.51*poly.getTimeMax());
	poly.solve();
    std::cout << "Total time to run polydisperse simulation: " << time.toc() << "s (Expected: 3s)" << std::endl;

    time.tic();
 	Contact highPoly(5.0);
 	highPoly.setTimeMax(0.44*highPoly.getTimeMax());
	highPoly.solve();
    std::cout << "Total time to run highly polydisperse simulation: " << time.toc() << "s (Expected: 3s)" << std::endl;

	return 0;
}

// create gperftools profile:

//#!/bin/bash
//		set +x
//		echo './profile.sh executable'
//
//# gperftools
//#https://github.com/ethz-asl/programming_guidelines/wiki/Profiling-Code
//#https://gperftools.github.io/gperftools/cpuprofile.html
//
//#create CPU profile
//CPUPROFILE=/tmp/cpu LD_PRELOAD=/usr/lib/libprofiler.so.0.4.5 $*
//		google-pprof --pdf $1 /tmp/cpu > /tmp/profile.pdf
//		gvfs-open /tmp/profile.pdf

// create gprof profile (you have to make with -pg)

//#!/bin/bash
//		set +x
//		echo './profile.sh executable'
//
//#https://linoxide.com/tools/gprof-performance-analysis-programs/
//
//#create CPU profile
//$*
//		gprof $1 gmon.out > /tmp/gprofile






