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

#include "DPMBase.h"
#include "Particles/SphericalParticle.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticSpecies.h"
#include <cmath>
#include <iostream>
#include <iomanip>


///This tests if moving the wall works with CWall::move(Vec3D velocity,Vec3D dt).
///The wall is moved in normal and tangential direction and the results 
///are compared with a system where the particles are moved instead of 
///the walls, with the same relative velocities
class WallSpecies : public DPMBase {
public:

	void setupInitialConditions() override {
		setXMax(20);
		setYMax(20);
		setZMax(20);
        setSystemDimensions(3);
        setParticleDimensions(3);
        setGravity(Vec3D(0,0,0));
        setTimeStep(.0002);
		setTimeMax(15);
        setSaveCount(500);

        //I don't know how to write this in a neat way.
        LinearViscoelasticSpecies* species0 = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        LinearViscoelasticSpecies* species1 = speciesHandler.copyAndAddObject(species0);
        LinearViscoelasticMixedSpecies* species01 = speciesHandler.getMixedObject(species0, species1);
        species0->setDensity(6.0/constants::pi);
        species0->setCollisionTimeAndRestitutionCoefficient(0.01,0.25,2);
        species01->setCollisionTimeAndRestitutionCoefficient(0.01,0.5,2);

		//set particles
		SphericalParticle P0;
        P0.setSpecies(speciesHandler.getObject(0));
		P0.setRadius(0.5);
		P0.setPosition(Vec3D(5,10,10));
		P0.setVelocity(Vec3D(-1,0,-1));
		particleHandler.copyAndAddObject(P0);
        
		//set walls
		InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
		w0.set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
		wallHandler.copyAndAddObject(w0);
        
        w0.set(Vec3D(-1,0,0), Vec3D(getZMin(), 0, 0));
        w0.setSpecies(speciesHandler.getObject(1));
		wallHandler.copyAndAddObject(w0);
        
 	}

    //Write the speed of the particles as the info line.
    double getInfo(const BaseParticle& P0) const override {
        return sqrt(mathsFunc::square(P0.getVelocity().X)+mathsFunc::square(P0.getVelocity().Z));
    }
};



int main(int argc UNUSED, char *argv[] UNUSED)
{
	std::cout << std::endl << "Simulation of one particle, d=0.1, interacting with walls with different coefficient of restitution" << std::endl;
	std::cout << std::endl << "The purpose of the test is to check walls species information is correctly picked up" << std::endl;
	WallSpecies problem;
	problem.setFileType(FileType::NO_FILE); //comment if you want file output
	problem.setName("WallSpecies");
	problem.solve();
	const Vec3D v = problem.particleHandler.getObject(0)->getVelocity();
	logger(INFO,"v_x(t_max)=% should be ~0.5",v.X);
	logger(INFO,"v_z(t_max)=% should be ~0.25",v.Z);
	helpers::check(v.Z,0.23642946655364,1e-7,"v_z");
	helpers::check(v.X,0.49298992225767,1e-7,"v_x");
}
