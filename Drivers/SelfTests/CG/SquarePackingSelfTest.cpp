//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#include <cmath>
#include <iostream>
#include <iomanip>
#include <Species/LinearViscoelasticSpecies.h>

class SquarePacking : public StatisticsVector<O> {

public:

	void setupInitialConditions() override {
        // make sure the number of particles is right
        logger(INFO, "Creating a cubic packing of %^3 particles", N);
        double Radius = .5;
        
        //set Particles' position, radius, velocity and bounding box
        setDomain({0, 0, -0.5}, {5, 5, 0.5});
        
        SphericalParticle P0;
        P0.setSpecies(speciesHandler.getObject(0));
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                
                P0.setRadius(Radius);
                P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
                P0.setPosition(Radius * Vec3D(1. + 2. * i, 1. + 2. * j, 0.0));
                //setXMax(std::max(getXMax(),P0.getPosition().X+Radius));
                //setYMax(std::max(getYMax(),P0.getPosition().Y+Radius));
                particleHandler.copyAndAddObject(P0);
                
            }	

		//set walls
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
		w0.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 1, 0, 0), Vec3D(getXMax(), 0, 0));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0,-1, 0), Vec3D(0, getYMin(), 0));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D( 0, 1, 0), Vec3D(0, getYMax(), 0));
		wallHandler.copyAndAddObject(w0);
	}
	
	int N;
	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    SquarePacking problem;
    problem.setName("SquarePackingSelfTest");
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    //set the number of particles
    problem.N = 5;
    problem.setSystemDimensions(2);
    problem.setParticleDimensions(2);
    species->setDensity(6.0 / constants::pi * 2. / 3.);
    //species->setDensity(1.9098593*2/3);
    problem.setGravity(Vec3D(0., -1, 0.));
    species->setCollisionTimeAndRestitutionCoefficient(.01, .1, 1.);
    problem.setTimeStep(.0002);
    problem.setTimeMax(1.0);
    problem.setSaveCount(1000);
    
    problem.setCGTimeMin(problem.getTimeMax() * .98);
    problem.setTimeMaxStat(1e20);
    
    logger(INFO, "Relax the packing");
    problem.solve();
    
    logger(INFO, "Get statistics");
    StatisticsVector<Y> stats("SquarePackingSelfTest");
    stats.statFile.setName("SquarePackingSelfTest.Y.stat");
    double n = 500;
    stats.setN(n);
    stats.setCGWidth(.1);
    stats.setCGTimeMin(problem.getTimeMax() * .98);
    stats.setTimeMaxStat(1e20);
    //stats.verbose();
    //stats.set_boundedDomain(true);
    stats.statistics_from_fstat_and_data();
    
    logger(INFO, "Get fully averaged statistics");
    StatisticsVector<O> statsO("SquarePackingSelfTest");
    statsO.statFile.setName("SquarePackingSelfTest.O.stat");
    statsO.setCGTimeMin(problem.getTimeMax() * .98);
    statsO.setTimeMaxStat(1e20);
    statsO.statistics_from_fstat_and_data();
    // should give you Density 1
}

