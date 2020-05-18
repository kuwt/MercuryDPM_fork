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

#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>
using constants::pi;

class ParticleInclusion : public Mercury3D
{
public:

	ParticleInclusion() {
		//name
		setName("ParticleInclusion");
        setMax(Vec3D(5,5,5));
        setMin(Vec3D(-5,-5,-5));
        setGravity(Vec3D(0,0,0));

        //contact properties
		LinearViscoelasticSpecies species;
		species.setDensity(6./pi);
		species.setStiffness(1);
		species.setDissipation(0);
        LinearViscoelasticSpecies* s = speciesHandler.copyAndAddObject(species);
        logger(INFO,"tc %", s->getCollisionTime(8));

        //add particles
        SphericalParticle p;
		p.setSpecies(s);
        p.setRadius(5);
        particleHandler.copyAndAddObject(p);

        p.setPosition(Vec3D(0.5,0,0));
        p.setRadius(1);
        particleHandler.copyAndAddObject(p);

        setTimeStep(0.001);
        setTimeMax(4.5);


        //xij'' = (ai-aj) = fij/mi + fij/mj = fij/mij = k/mij*(6-xij)
        //mij = 125*8/(125+8) approx 8
        //x'' = 3/4 - x/8
	}

    void printTime() const override
    {
        logger(INFO,"t % EKin % EEla % C %",getTime(),getKineticEnergy(),getElasticEnergy(),interactionHandler.getLastObject()->getContactPoint());
    }


};

int main(int argc, char *argv[])
{
	ParticleInclusion ParticleInclusion;
	ParticleInclusion.solve(argc,argv);
}
