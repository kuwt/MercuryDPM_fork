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

//#include "Sinter2.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <Species/LinearPlasticViscoelasticSlidingFrictionSpecies.h>
#include <Particles/BaseParticle.h>
#include <Chute.h>

class two_particle_collision : public Chute {
public:
	two_particle_collision () : Chute() {
        //set default species
		Mdouble d = 4.5e-6;
        species->setDensity(1000);
		Mdouble m = 4.0/3.0*constants::pi*species->getDensity()*d*d*d;
		Mdouble tc = constants::pi/sqrt(1344/(m/2.0));
		Mdouble r = 0.88;
		Mdouble g = 9.8e3;
		N=1;
		M=1;

        species->setSlidingFrictionCoefficient(0);
	    Mdouble k2_k1_ratio = 5;
	    Mdouble kc_k1_ratio = 5;
	    Mdouble depth = 0.05; //delta0max=d/10=4.5e-7
		species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r, r, m);
	    setInflowParticleRadius(d/2.0);
	    setFixedParticleRadius(d/2.0);
        species->setPlasticParameters(species->getLoadingStiffness()/k2_k1_ratio, species->getLoadingStiffness(), species->getLoadingStiffness()/k2_k1_ratio*kc_k1_ratio, depth);
        setTimeStep(tc / 20.0);
		setTimeMax(tc*10e4);
	    setSaveCount(getTimeMax()/ getTimeStep()/100.0);
	    setGravity(Vec3D(0,0,-g));
		setName("Bed");

        speciesHandler.copyAndAddObject(speciesHandler.getObject(0));
		species->setPenetrationDepthMax(2.0*species->getPenetrationDepthMax());
        species->setLoadingStiffness(0.5 * species->getLoadingStiffness());
        
        nCreated_=0;
	}

	void setupInitialConditions() override
	{
		setXMax(N*2.0*getInflowParticleRadius());
		setYMax(N*2.0*getInflowParticleRadius());
		setZMax(M*2.0*getInflowParticleRadius());
		
		wallHandler.clear();
		InfiniteWall w;
		w.set(Vec3D(0, 0, -1), Vec3D(0,0,getZMin()-getFixedParticleRadius()));
		wallHandler.copyAndAddObject(w);
		
		boundaryHandler.clear();
		PeriodicBoundary b;
		b.set(Vec3D(1, 0, 0), getXMin(), getXMax());
		boundaryHandler.copyAndAddObject(b);
		b.set(Vec3D(0, 1, 0), getYMin(), getYMax());
		boundaryHandler.copyAndAddObject(b);

		createBottom();
		add_flow_particles();
		writeRestartFile();
		write(std::cout,false);
	}
	
		//add flow particles
	void add_flow_particles() 
	{
		Mdouble Nbase = particleHandler.getNumberOfObjects();
		particleHandler.setStorageCapacity(N*N*M+Nbase);
		//setHGridNumberOfBucketsToPower(N*N*M);
		hGridActionsBeforeTimeLoop();
		hGridActionsBeforeTimeStep();
		setInflowHeight(2.0*getZMax());

		writeRestartFile();
		//try to find new insertable particles
		while (particleHandler.getNumberOfObjects()<N*N*M+Nbase){
			create_inflow_particle();
            if (checkParticleForInteraction(P0))
            {
                particleHandler.copyAndAddObject(P0);
                increaseNCreated();
            }
            else
            {
                setInflowHeight(getInflowHeight() + .0001* getMaxInflowParticleRadius());
            }

		}
		//setHGridNumberOfBucketsToPower();
	}
	
	//defines type of flow particles	
	void create_inflow_particle()
	{
		P0.setSpecies(speciesHandler.getObject(1));
		P0.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
		//P0.computeMass();
		
		
        //The position components are first stored in a Vec3D, because if you pass them directly into setPosition the compiler is allowed to change the order in which the numbers are generated
        Vec3D position;
		position.X = random.getRandomNumber(getXMin()+P0.getRadius(),getXMax()-P0.getRadius());
		position.Y = random.getRandomNumber(getYMin()+P0.getRadius(),getYMax()-P0.getRadius());
		position.Z = random.getRandomNumber(getZMin()+2.0*P0.getRadius(), getInflowHeight());
        P0.setPosition(position);
		P0.setVelocity(Vec3D(0.0,0.0,0.0));
	}

	// void actionsBeforeTimeStep() override {
	// 	//std::cout << getTime() << std::endl;
	// 	if (getTime()<0.3*getTimeMax())
	// 		set_Temperature(0);
	// 	else if (getTime()<0.7*getTimeMax())
	// 		set_Temperature(1.3);
	// 	else
	// 		set_Temperature(0);
	// }		
	
    int getNCreated() const
    {
        return nCreated_;
    }

    void increaseNCreated()
    {
        nCreated_++;
    }

    int nCreated_;
    //causes segmentation faults if this is not there;
    SphericalParticle P0;
public:
	int N, M;
    LinearPlasticViscoelasticSlidingFrictionSpecies* species;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	two_particle_collision md;
    md.N=10;
    md.M=5;
	md.setName("Bed");
	Mdouble R = md.getInflowParticleRadius();
	md.setInflowParticleRadius(0.96*R,R);
    md.setRoughBottomType(MONOLAYER_DISORDERED);
	md.setXBallsAdditionalArguments("-v0 -solidf");
 	md.fStatFile.setFileType(FileType::NO_FILE);
 	md.setHGridMaxLevels(1);
    if (argc>=2)
		md.M = atoi(argv[1]);
	else
		md.M = 5;
	std::cout << "H" << md.M << " W" << md.N << std::endl;
	std::stringstream s;
	s << "BedH" << md.M;
	md.setName(s.str().c_str());
	md.solve();
}
//to check if the bed is relaxed:
//set logscale y
//gnuplot> p 'Bed.ene' u 1:($3/$5)
