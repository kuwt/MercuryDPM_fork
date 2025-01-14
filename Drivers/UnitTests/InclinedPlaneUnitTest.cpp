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

#include "DPMBase.h"
#include "Particles/SphericalParticle.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>

class InclinedPlane : public DPMBase {
public:
	
	void setupInitialConditions() override {
		InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
		w0.set(Vec3D(0,-1,0),Vec3D(0, getYMin(), 0));
		wallHandler.copyAndAddObject(w0);
		
		SphericalParticle p0; //important
        p0.setSpecies(speciesHandler.getObject(0));
		p0.setRadius(0.005);
		p0.setPosition(Vec3D((getXMin()+getXMax())/2.0,getYMin()+p0.getRadius(),0));
		p0.setVelocity(Vec3D(0.0,0.0,0.0));
		particleHandler.copyAndAddObject(p0);
	}

};

int main(int argc, char *argv[])
{
	logger(INFO, "Single particle rolling on a bottom plate");
	///Start off my solving the default problem
	InclinedPlane problem;
    
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    problem.setName("InclinedPlane");
    problem.setFileType(FileType::NO_FILE); //comment if you want file output
    problem.setSaveCount(50000);
    problem.setSystemDimensions(2);
	problem.setParticleDimensions(3);
	Mdouble theta = 15.0*constants::pi/180.0;
	problem.setGravity(9.8*Vec3D(std::sin(theta),-std::cos(theta),0.0));
	species->setDensity(2000);
    species->setStiffness(10000);
    species->setDissipation(0);
	species->setSlidingStiffness(10000);
	species->setSlidingDissipation(0);
	species->setSlidingFrictionCoefficient(0.5);
	//species->setRollingStiffness(10000);
	//species->setRollingDissipation(0);
	//species->setRollingFrictionCoefficient(0.5);
	//tc=7e-5
    problem.setTimeStep(0.02*species->getCollisionTime(species->getMassFromRadius(0.005)));
	problem.setTimeMax(1.0);
	problem.setXMax(0.1);
	problem.setYMax(0.1);
	problem.setZMax(0.1);
	problem.solve(argc, argv);
	//problem.write(std::cout, true);
    const Vec3D v = problem.particleHandler.getObject(0)->getVelocity();
    logger(INFO,"v_x(t_max)=% should be ~1.81",v.X);
    logger(INFO,"v_y(t_max)=% should be ~0",v.Y);
    helpers::check(v.X,1.81167,1e-4,"v_x");
    helpers::check(v.Y,0.00234107,1e-4,"v_y");

}
