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

//choose rvar (\delta_a^T) small, f.e. 1e-4 K^{-1} 
#include <iomanip>

#include "Chute.h"
#include "Particles/BaseParticle.h"
#include <iostream>
#include "Species/LinearPlasticViscoelasticSlidingFrictionSpecies.h"

class Sinter : public Chute{
public:

	void add_particles() {}

	Sinter () : Chute() {
        Temperature = 0;	
		NeckGrowthSpeed = 1e-5;	//neck grows like (1\pm\exp(-k*t))*FinalNeckSize
		
		//set default species
		Mdouble m = 1;
		Mdouble d = 1;
		Mdouble g = 1;
		Mdouble tc = 0.05;
		Mdouble r = 0.88;
        species->setSlidingFrictionCoefficient(0);
	    Mdouble k2_k1_ratio = 5;
	    Mdouble kc_k1_ratio = 5;
	    Mdouble depth = 0.2;
	    
		species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r, r, m);
	    setInflowParticleRadius(d/2.0);
	    setFixedParticleRadius(d/2.0);
        species->setDensity(6./constants::pi*m/d/d/d);
	    species->setPlasticParameters(species->getLoadingStiffness()/k2_k1_ratio, species->getLoadingStiffness(), species->getLoadingStiffness()/k2_k1_ratio*kc_k1_ratio, depth);
	    setTimeStep(tc/50.0);
		setTimeMax(tc*4000.0);
	    setSaveCount(400);
	    setGravity(Vec3D(0,0,-g));
		setName("Sinter");
        *oldSpecies=*species;

    }

	void setupInitialConditions()
	{
        *oldSpecies=*species;
    }

 	void printTime() const {
	    std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
            << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
            << ", T=" << std::setprecision(3) << std::left << std::setw(6) << Temperature
            << ", k1=" << std::setprecision(3) << std::left << std::setw(6) << species->getLoadingStiffness()
            << "\n";
		std::cout.flush();
	}

  	void setTemperature(Mdouble)
    { 
        for (auto it=interactionHandler.begin(); it!=interactionHandler.end(); it++)
            if ((*it)->getTimeStamp()>=getTime())
            {
                auto interaction = dynamic_cast<Interaction<LinearPlasticViscoelasticInteraction,SlidingFrictionInteraction>*>(*it);
                interaction->setMaxOverlap(interaction->getMaxOverlap()+NeckGrowthSpeed*getTimeStep());
            }
	}

	///Interaction parameters of the sintered material in cold state
	Mdouble Temperature; //determines final neck size (up to k2max limit)
	Mdouble NeckGrowthSpeed; //in [m/s]
    LinearPlasticViscoelasticSlidingFrictionSpecies* species;
    LinearPlasticViscoelasticSlidingFrictionSpecies* oldSpecies;
};
