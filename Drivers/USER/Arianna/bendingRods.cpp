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


#include <iostream>
#include "Species/LinearViscoelasticFrictionChargedBondedSpecies.h"
#include "Mercury3D.h"
using constants::pi;

/**
 * Tests the bending behaviour of bonded particles.
 * A chain of n particles is created, the first particle is fixed, then the chain is bent down by gravity.
 * Finally, at t=0.5t_{Max}, gravity is zet to zero, the the chain relaxes again, almost retrieving its original position.
 */
class BendingRods : public Mercury3D {
public:
    // Constructor, takes the number of particles in the chain as input
    BendingRods (unsigned n)
    {
        logger.assert_always(n>1,"The number of particles in the chain has to be at least two");

        //set name, gravity, output options
        setName("bendingRods");
        setGravity({0,0,-5e4});
        //setXBallsAdditionalArguments("-solidf -v0 -cmode 8");
        setTimeMax(1.0);

        //Create the same contact properties as in RodsEF_2D.cpp (except dissipation)
        LinearViscoelasticFrictionChargedBondedSpecies s;
        s.setDensity(6./pi);
        s.setStiffness(1e6);
        s.setSlidingStiffness(2./7.*s.getStiffness());
        s.setSlidingDissipation(2./7.*s.getDissipation());
        s.setSlidingFrictionCoefficient(std::numeric_limits<Mdouble>::infinity());
        s.setRollingStiffness(2./5.*s.getStiffness());
        s.setRollingDissipation(2./5.*s.getDissipation());
        s.setRollingFrictionCoefficient(std::numeric_limits<Mdouble>::infinity());
        s.setBondForceMax(4e5);
        auto species = speciesHandler.copyAndAddObject(s);

        // create a chain of particles
        BaseParticle p;
        p.setSpecies(species);
        p.setRadius(0.5);
        for (unsigned i = 0; i<n; i++)
        {
            particleHandler.copyAndAddObject(p);
            p.setPosition(p.getPosition()+Vec3D(1.2*p.getRadius(),0,0));
        }
        // fix the first particle
        particleHandler.getObject(0)->fixParticle();
        // bond the chain
        for (unsigned i = 1; i<n; i++)
        {
            dynamic_cast<ChargedBondedInteraction *>(
             interactionHandler.getInteraction(particleHandler.getObject(i-1),particleHandler.getObject(i), 0)
            )->bond();
        }

        //set the time step according to the species properties
        Mdouble mass = species->getMassFromRadius(p.getRadius());
        Mdouble tc = species->getCollisionTime(mass);
        Mdouble r = species->getRestitutionCoefficient(mass);
        logger(INFO,"Collision time %",tc);
        logger(INFO,"Restitution coefficient %",r);
        setTimeStep(0.02*tc);
        setSaveCount(1000);
        eneFile.setSaveCount(50);
        restartFile.setSaveCount(std::numeric_limits<int>::max());
        fStatFile.setFileType(FileType::NO_FILE);
        setTimeMax(2000*tc);

        //fix the domain size to the size of the chain
        setMin(p.getRadius()*Vec3D(-1,-1,-1));
        setMax(p.getRadius()*Vec3D(1,1,1)+particleHandler.getLastObject()->getPosition());
    }

    // add background dissipation (simpler, as is dissipates energy quickly, but unrealistic)
    void computeExternalForces(BaseParticle* CI) override
    {
        DPMBase::computeExternalForces(CI);
        CI->addForce(-650*CI->getVelocity());
        CI->addTorque(-0.1*650*CI->getAngularVelocity()*CI->getRadius());
    }

    // change gravity after half the time
    void actionsBeforeTimeStep() override
    {
        if (getTime()<0.5*getTimeMax())
        {}//setGravity({0,0,-45e4*sin(2.*pi*getTime()/getTimeMax())});
        else
            setGravity({0,0,0});

    }

    // write some outout to the ene file to test the behaviour of the tangential springs
    // For example, check that the rolling spring is orthogonal to the normal contact direction by plotting in gnuplot:
    // p 'bendingRods.ene' u ($5*$8+$7*$10) w lp
    void writeEneHeader(std::ostream& os) const override
    {
        os << "Time\tSlidingSpring\t\tRollingSpring\t\tPosition\t\tAngularVelocity\n";
    }

    // write some outout to the ene file to test the behaviour of the tangential springs
    void writeEneTimeStep(std::ostream& os) const override
    {
        auto i = dynamic_cast<const FrictionInteraction*>(interactionHandler.getObject(0));
        auto p = particleHandler.getObject(1);
        os << getTime() << '\t' << i->getSlidingSpring() << '\t' << i->getRollingSpring() << '\t'
           << p->getPosition() << '\t' << p->getAngularVelocity() << '\n';
    }

};

int main () {
    BendingRods br(5);
    br.solve();
    return 0;
}

























