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
#include <Species/LinearViscoelasticSlidingFrictionBondedSpecies.h>
#include <cassert>
class SlidingFrictionUnitTest : public DPMBase {

	void writeEneHeader(std::ostream& os) const override
	{
		os << "time slidingSpringX slidingSpringY slidingSpringZ\n";
	}

	void writeEneTimeStep(std::ostream& os) const override
	{
		auto i = dynamic_cast<const LinearViscoelasticSlidingFrictionBondedSpecies::InteractionType*>(interactionHandler.getLastObject());
		assert(i);
		os << getTime() << ' ' << i->getSlidingSpring() << '\n';
	}
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	//Unit diameter
	Mdouble radius = 0.5;
	//Equilibrium overlap
	Mdouble overlap = 0.1*radius;
	//Collision time (determines time step)
	Mdouble tc = 50e-3;
	//rotation velocity
	Mdouble velocity = 1.0;

	SlidingFrictionUnitTest sf;
	sf.setGravity({0,0,0});
	sf.setTimeStep(0.02*tc);
	sf.setName("SlidingFrictionUnitTest");
	sf.setMin(-radius*Vec3D(2,2,2));
	sf.setMax( radius*Vec3D(2,2,2));
	sf.setDimension(3);
	sf.setTimeMax(2.0);
	sf.eneFile.setSaveCount(1);

	//Add species
	auto s = sf.speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionBondedSpecies());
	//Set unit density
	s->setDensity(1.0);
	//Make particles relatively soft and dissipative
	Mdouble mass = s->getMassFromRadius(radius);
	s->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc,0.5,1.0,mass);
	//Set a high friction, so sliding movement is mostly elastic
	s->setSlidingFrictionCoefficient(1e20);
	//10% overlap in equilibrium position
	s->setBondForceMax(overlap*s->getStiffness());
	// check that stiffness is not too small
	logger(INFO,"Stiffness k=%",s->getStiffness());

	//Add particles
	SphericalParticle p;
	p.setSpecies(s);
	p.setRadius(radius);
	p.setPosition({ radius-0.5*overlap,0,0});
	p.setVelocity(velocity*Vec3D(0,0,1));
	auto p0 = sf.particleHandler.copyAndAddObject(p);
	p.setPosition({-radius+0.5*overlap,0,0});
	p.setVelocity(velocity*Vec3D(0,0,-1));
	auto p1 = sf.particleHandler.copyAndAddObject(p);

	//Add bonded interaction
	auto i = dynamic_cast<LinearViscoelasticSlidingFrictionBondedSpecies::InteractionType*>(sf.interactionHandler.getInteraction(p0,p1,0));
	i->bond();
	i->setSlidingSpring({0,0,overlap});
	sf.solve();

	helpers::writeToFile("SlidingFrictionUnitTest.gnu","set size ratio -1; p 'SlidingFrictionUnitTest.ene' u 2:4 w lp");
}
