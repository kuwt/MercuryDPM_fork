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

//choose rvar (\delta_a^T) small, f.e. 1e-4 K^{-1} 
#include <iomanip>
#include <iostream>

#include "Chute.h"
#include "Particles/BaseParticle.h"
#include "Walls/BaseWall.h"
#include "Species/Species.h"
#include "Species/FrictionForceSpecies/SlidingFrictionSpecies.h"
#include "LinearSinterSpecies.h"

class Sinter : public Chute{
public:

	void add_particles() {}

	Sinter () : Chute()
    {
        species = speciesHandler.copyAndAddObject(Species<LinearSinterSpecies,SlidingFrictionSpecies>());
        setName("Sinter");
        setSaveCount(400);
        setGravity(Vec3D(0,0,0));
        scaleSpeciesProperties(1.0, 1.0);
    }

    void scaleSpeciesProperties(Mdouble diameter, Mdouble density)
    {
		//set default species
		Mdouble mass = constants::pi/6.0*density*mathsFunc::cubic(diameter);
		Mdouble tc = 0.05;
        Mdouble eps = 0.88;
        Mdouble beta = eps;
        Mdouble k2_k1_ratio = 5.0;
	    Mdouble kc_k1_ratio = 5.0;
	    Mdouble depth = 0.2;

        //species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, r, r, m);
        species->setDissipation(-mass / tc * std::log(eps));
        species->setLoadingStiffness(.5 * mass * (mathsFunc::square(constants::pi/tc) + mathsFunc::square(species->getDissipation()) /mass));
        species->setSlidingStiffness(2.0 / 7.0 * species->getLoadingStiffness() * (mathsFunc::square(constants::pi) + mathsFunc::square(log(beta))) / (mathsFunc::square(constants::pi) + mathsFunc::square(log(eps))));
        if (beta != 0.0)
            species->setSlidingDissipation(-2 * log(beta) * sqrt(1.0 / 7.0 * mass * species->getSlidingStiffness() / (mathsFunc::square(constants::pi) + mathsFunc::square(log(beta)))));
        else
            species->setSlidingDissipation(2. * sqrt(1.0 / 7.0 * mass * species->getSlidingStiffness()));
        species->setPlasticParameters(species->getLoadingStiffness()/k2_k1_ratio, species->getLoadingStiffness(), species->getLoadingStiffness()/k2_k1_ratio*kc_k1_ratio, depth);

        maxNeckGrowthRate = depth*diameter/tc; // maximum growth after tc/Temperature;
        species->setDensity(density);
        setInflowParticleRadius(diameter/2.0);
        setFixedParticleRadius(diameter/2.0);
        
        setTimeStep(tc/50.0);
		setTimeMax(tc*4000.0);
    }

    Species<LinearSinterSpecies,SlidingFrictionSpecies>* species;
    Mdouble maxNeckGrowthRate;
};
