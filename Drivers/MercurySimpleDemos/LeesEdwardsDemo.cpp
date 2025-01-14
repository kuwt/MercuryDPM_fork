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

#include "Mercury2D.h"
#include "Boundaries/LeesEdwardsBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Particles/BaseParticle.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"

class LeesEdwardsDemo : public Mercury2D
{
    public:

        LeesEdwardsDemo()
        {
            setName("LeesEdwardsDemo");
            dataFile.setFileType(FileType::MULTIPLE_FILES);
            fStatFile.setFileType(FileType::NO_FILE);

            r = 0.05;
            Mdouble mass = M_PI * r * r;

            auto species = speciesHandler.copyAndAddObject(new LinearViscoelasticFrictionSpecies);
            species->setDensity(1);
            species->setCollisionTimeAndRestitutionCoefficient(2e-2, 0.3, mass);
            logger(WARN, "collision time %, restitution coefficient %, max vel %",
                    species->getCollisionTime(mass), species->getRestitutionCoefficient(mass),
                    species->getMaximumVelocity(r, mass));
            species->setSlidingFrictionCoefficient(tan( 21.8 * M_PI / 180. ));
            species->setSlidingStiffness(2.0/7.0 * species->getStiffness());
            species->setSlidingDissipation(2.0/7.0 * species->getDissipation());
            species->setRollingFrictionCoefficient(tan( 0.0 * M_PI / 180. ));
            species->setRollingStiffness(2.0/5.0 * species->getStiffness());
            species->setRollingDissipation(2.0/5.0 * species->getDissipation());
            species = speciesHandler.copyAndAddObject(species);

            //set time and file properties
            setTimeStep(species->getCollisionTime(mass) / 20.0);
            setTimeMax(25.0);
            setSaveCount(0.2 / getTimeStep());

            setMin({0, -1, 0});
            setMax({3, 1, 1});
            setGravity({0, 0, 0});
        }

        void setupInitialConditions() override
        {
            auto species = speciesHandler.getObject(0);

            const Mdouble v = 2.0;
            LeesEdwardsBoundary leesEdwardsBoundary;
            leesEdwardsBoundary.set(
                    [v](double t) { return v*t; },
                    [v](double t) { return v; },
                    getXMin(), getXMax(), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(leesEdwardsBoundary);

            //define common particle properties
            auto p = new SphericalParticle;
            p->setSpecies(species);
            const Mdouble volumeFraction = 0.7;
            logger(INFO,"Inserting particles of diameter % to %, volumeFraction %",
                    .9*r, 1.1*r, volumeFraction);
            auto insb = boundaryHandler.copyAndAddObject(new CubeInsertionBoundary);
            Mdouble vel = 0.1;
            insb->set(p, 1,
                    Vec3D(getXMin(), getYMin(), 0), Vec3D(getXMax(), getYMax(), 0),
                    Vec3D(-vel, -vel, 0), Vec3D(vel, vel, 0), 0.9*r, 1.1*r);
            insb->checkBoundaryBeforeTimeStep(this);
        }

        void actionsAfterTimeStep() override
        {
            if (boundaryHandler.getNumberOfObjects() == 2)
            {
                Mdouble volfrac = particleHandler.getVolume() / (
                        (getXMax() - getXMin()) * (getYMax() - getYMin())) ;
                if (volfrac > 0.7)
                {
                    logger(DEBUG, "Finished inserting particles");
                    boundaryHandler.removeObject(1);
                    for (auto p : particleHandler)
                        p->setInfo(p->getPosition().Y);
                }
                else
                    logger(DEBUG, "volumefraction = %", volfrac);
            }
        }

        Mdouble r;
};

int main()
{
    //instantiate the class
    LeesEdwardsDemo problem;
    //set output and time stepping properties
    //problem.setXBallsAdditionalArguments("-w0 -v0 -solidf -cmode 5");
    //solve
    problem.solve();
    return 0;
}

