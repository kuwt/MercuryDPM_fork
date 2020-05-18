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
#include <Species/LinearViscoelasticFrictionSpecies.h>
using constants::pi;


/**
 * This code tests the torsion force model, as published in Luding 2008.
 * A particle is squeezed in between two fixed particles, and rotated such that a torsional torque is created.
 * You should see that the torque is capped at torsionfriction*radius*normalForce, then oscillates and decays.
 */
class DPM : public Mercury3D{
public:

    DPM()
    {
        //set parameters
        Mdouble collisionTime = 2e-2;
        Mdouble restitution = 0.7;
        Mdouble torsionFriction = 0.1;
        Mdouble N = 1000;

        //define species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        species->setDensity(6.0 / pi);
        species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(collisionTime,restitution,restitution,1);
        species->setTorsionStiffness(0.4*species->getStiffness());
        species->setTorsionDissipation(0.4*species->getDissipation());
        species->setTorsionFrictionCoefficient(torsionFriction);

        setTimeStep(0.02*collisionTime);
        setTimeMax(10*collisionTime);
        setSaveCount(5);

        setName("TorsionFrictionSelfTest");
        setXBallsAdditionalArguments("-v0 -solidf -3dturn 1");

        SphericalParticle p(species);
        p.setRadius(0.5);

        p.setPosition({0,0,0});
        p.setAngularVelocity({0,0,1});
        particleHandler.copyAndAddObject(p);

        p.fixParticle();
        p.setAngularVelocity({0,0,0});

        p.setPosition({0,0,-1.99*p.getRadius()});
        particleHandler.copyAndAddObject(p);
        p.setPosition({0,0,1.99*p.getRadius()});
        particleHandler.copyAndAddObject(p);

        setDomain({-1,-1,-1},{1,1,1});
    }

    void actionsAfterTimeStep() override
    {
        static std::ofstream ofile ("TorsionFrictionSelfTest.torque");
        auto p = particleHandler.getLastObject();
        ofile << getTime()
              << ' ' << p->getTorque().Z
              << ' ' << p->getForce().Z
              << '\n';
    }

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    DPM dpm;
    dpm.solve();

    logger(INFO,"Run gnuplot %.gnu --persist to view output; you should see that the torque is capped at torsionfriction*radius*normalForce, then oscillates and decays",dpm.getName());
    helpers::writeToFile(dpm.getName()+".gnu",
                         "set xlabel 'time'\n"
                         "set ylabel 'torque'\n"
                         "p 'TorsionFrictionSelfTest.torque' u 1:2 w l t 'torque', '' u 1:(0.1*0.5*$3) w l t 'torsionfriction*radius*normalForce'"
    );
}
