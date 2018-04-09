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

#include "Sinter.h"
#include <iostream>
#include "Particles/BaseParticle.h"

class two_particle_collision : public Sinter {
public:
	two_particle_collision () : Sinter() {
		//set default species
        setName("SinterParticlePair");
        //setSaveCount(1);
        Mdouble diameter = 4.5e-6;
        Mdouble density = 1000.0;
		scaleSpeciesProperties(diameter, density);
	}

	void setupInitialConditions()
	{
		//calculate a deltaMax and a corresponding delta0
 		Mdouble R = getInflowParticleRadius();

        Mdouble delta0Star = species->getPenetrationDepthMax() *2.0*R;
        Mdouble deltaMaxStar = (species->getUnloadingStiffnessMax() /(species->getUnloadingStiffnessMax() -species->getLoadingStiffness())) * delta0Star;
        Mdouble deltaMax = 0.0005*R;
        Mdouble l = deltaMax/deltaMaxStar;
        Mdouble k2 = species->getLoadingStiffness() +l*(species->getUnloadingStiffnessMax() -species->getLoadingStiffness());
        Mdouble delta0 = (k2-species->getLoadingStiffness())/k2*deltaMax;


        //create particles with overlap delta0
		BaseParticle P;
		P.setRadius(R);
		P.setPosition(Vec3D(R,R,-R+delta0/2.0));
        BaseParticle* PI = particleHandler.copyAndAddObject(P);

        P.setPosition(Vec3D(R,R, R-delta0/2.0));
        BaseParticle* PJ = particleHandler.copyAndAddObject(P);

        BaseInteraction* C = PJ->getInteractionWith(PI, getTime(), &interactionHandler);
        dynamic_cast<LinearSinterInteraction*>(C)->setMaxOverlap(deltaMax);
        std::cout << "C" << interactionHandler.getNumberOfObjects() << std::endl;
        species->setNeckGrowthRate(0);
        std::cout << "C" << species->getNeckGrowthRate() << std::endl;

        setXMax(2.0*R);
		setYMax(2.0*R);
		setZMin(-2.0*R);
        setZMax(2.0*R);

        writeRestartFile();
		write(std::cout,true);
	}
	
	void actionsBeforeTimeStep(){
        if (getTime()+ getTimeStep()>0.1*getTimeMax() && getTime()<=0.1*getTimeMax())
            species->setNeckGrowthRate(0.0001*maxNeckGrowthRate);
        if (getTime()+ getTimeStep()>0.9*getTimeMax() && getTime()<=0.9*getTimeMax())
            species->setNeckGrowthRate(0);
	}
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	///Start off by solving the default problem
	two_particle_collision md;
	//md.setupInitialConditions();
    md.solve();
}
//now plot the position of one particle:
//p 'SinterParticlePair.data' u ($3*1e6) every 3::1
