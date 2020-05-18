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
#include "Logger.h"

/*!
 * Simple contacts, without any walls, boundaries, output, hGrid
 */
class Wall : public Mercury3D
{
	Mdouble polydispersity_;

public:

    explicit Wall ()
	{
		LinearViscoelasticSpecies s;
		s.setDensity(2000);
		s.setStiffness(1e3);
		speciesHandler.copyAndAddObject(s);

		setTimeStep(5e-5);
	    setTimeMax(300);
    	setName("SpeedTestWallInteractions");
		setFileType(FileType::NO_FILE);
		//dataFile.setFileType(FileType::ONE_FILE);
	}

	void setupInitialConditions() override
	{
        Mdouble r = 0.002;
        Mdouble h = 1.0*r;
		setDomain({-h,-h,-h},{h,h,h});
    	setSaveCount(NEVER);

        const LinearViscoelasticSpecies* s = dynamic_cast<const LinearViscoelasticSpecies*>(speciesHandler.getObject(0));
		SphericalParticle p;
        p.setSpecies(s);
		p.setRadius(r);
        double v = 0.1;
		p.setVelocity({v,0,0}); //expected overlap: mv^2=kdel^2 -> del = sqrt(m/k) v
        double k = s->getStiffness();
        double m = p.getMass();
        double tc = s->getCollisionTime(m);
        double deltaMax = sqrt(m/k)*v;
        logger(INFO,"m % k % v % tc % deltaMax/r %",m,k,v,tc,deltaMax/r);
		particleHandler.copyAndAddObject(p);

		InfiniteWall w;
		w.setSpecies(s);
        w.set({1,0,0},{h,0,0});
        wallHandler.copyAndAddObject(w);
        w.set({-1,0,0},{-h,0,0});
        wallHandler.copyAndAddObject(w);
	}
};

int main()
{
	std::cout << "An elastic particle between two walls." << std::endl;
	
    Time time;
    Wall wall;
	wall.solve();
    std::cout << "Total run time: " << time.toc() << "s (Expected: 3s)" << std::endl;
	return 0;
}
