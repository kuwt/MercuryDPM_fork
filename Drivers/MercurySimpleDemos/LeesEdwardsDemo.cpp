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

#include "Mercury2D.h"
#include "Boundaries/LeesEdwardsBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"

class LeesEdwardsDemo : public Mercury2D
{
public:

    void setupInitialConditions()
    {
        dataFile.setFileType(FileType::MULTIPLE_FILES);

        //set parameters to define the species properties
        const Mdouble particleRadius=0.04; // such that diameter is one
        const Mdouble collisionTime = 6.3291e-3; //relatively stiff particles
        const Mdouble restitution = 0.1; //restitution close to glass particles
        const Mdouble sizeDistribution = 1.1;
        const Mdouble volumeFraction = 0.7;
        const Mdouble velocity = 1.0;

        //define species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(4.0/constants::pi); //such that particle mass is one
        species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(
         collisionTime, restitution, restitution, species->getMassFromRadius(particleRadius)); //set stiffness and dissipation
        species->setSlidingFrictionCoefficient(0.5);

        //set gravity, time step
        setGravity(Vec3D(0,0,0)); // such that gravity is one
            setTimeMax(1000);
            setTimeStep(1.6e-4);
            setSaveCount(625);

        //set domain size
        setXMin(0); setYMin(0); setZMin(0);
        setXMax(2); setYMax(1); setZMax(1);

        //define leesEdwardsBoundary
        LeesEdwardsBoundary leesEdwardsBoundary;
        leesEdwardsBoundary.set(
         [velocity](double time) { return time * velocity; },
         [velocity](double time UNUSED) { return velocity; },
         getXMin(), getXMax(), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(leesEdwardsBoundary);

        //define common particle properties
        auto p = new SphericalParticle;
        p->setSpecies(speciesHandler.getObject(0));
        p->setRadius(particleRadius);
        Mdouble rMin=2.0*particleRadius/(sizeDistribution+1);
        Mdouble rMax=sizeDistribution*rMin;
        p->setRadius(rMax);
        logger(INFO,"Inserting particles of diameter % to %, volumeFraction %", 2.0*rMin, 2.0*rMax, volumeFraction);

        Vec3D position = {0,0,0};

        PSD myPSD;
        myPSD.setDistributionUniform(rMax,rMax,10000);
        auto insb = boundaryHandler.copyAndAddObject(new CubeInsertionBoundary);
        Mdouble vel = 0.1;
        insb->set(p, 1, 
                Vec3D(getXMin(), getYMin(), 0), Vec3D(getXMax(), getYMax(), 0),
                Vec3D(-vel, -vel, 0), Vec3D(vel, vel, 0));
        insb->setPSD(myPSD);
        insb->checkBoundaryBeforeTimeStep(this);
    }

    void actionsAfterTimeStep()
    {
            if (boundaryHandler.getNumberOfObjects() == 2)
            {
                if (particleHandler.getVolume() > 0.7 * (getXMax() - getXMin()) * (getYMax() - getYMin()))
                {
                    logger(WARN, "Finished inserting particles");
                    boundaryHandler.removeObject(1);
                    for (auto p : particleHandler)
                        p->setInfo(p->getPosition().Y);
                }
                else
                    logger(WARN, "volumefraction = %", particleHandler.getVolume()/2);
            }
    }


};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    //instantiate the class
    LeesEdwardsDemo leesEdwards;
    //set name
    leesEdwards.setName("LeesEdwardsDemo");
    //set output and time stepping properties
    leesEdwards.setXBallsAdditionalArguments("-w0 -v0 -solidf -cmode 5");
    //solve
    leesEdwards.solve();
    return 0;
}

