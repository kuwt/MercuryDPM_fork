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

#include "Mercury3D.h"
#include "Particles/SuperQuadricParticle.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Math/ExtendedMath.h"

class SphericalSuperQuadricCollision : public Mercury3D
{
    void setupInitialConditions() override
    {
    
        LinearViscoelasticSpecies species;
        species.setStiffness(2e5);
        species.setDissipation(25);
        speciesHandler.copyAndAddObject(species);
        
        SuperQuadricParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        //temporary hack:
        p0.setAxesAndExponents(1.0,1.0,1.0,1.0,1.0);
        p0.setInertia();
        SuperQuadricParticle p1 = *(p0.copy());
    
        p0.setPosition(Vec3D(0, 0.0, 0.0));
        p1.setPosition(Vec3D(3, 0.0, 0.0));
        p0.setVelocity(Vec3D(1, 0.0, 0.0));
        p1.setVelocity(Vec3D(-1, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);
        particleHandler.copyAndAddObject(p1);
        logger(INFO, "interaction radius p0: %", p0.getMaxInteractionRadius());
    
    
        setTimeStep(species.getCollisionTime(1) / 50);
        setTimeMax(1);
        setMin(-1, -1, -1);
        setMax(4, 1, 1);
    }
    
    //check contact time
    void actionsAfterTimeStep() override
    {
        BaseInteraction* interaction = particleHandler.getObject(0)->getInteractionWith(particleHandler.getObject(1), getTimeStep(), &interactionHandler);
        if (interaction!=nullptr)
        {
            logger.assert_always(getTime() > 0.499 && getTime() < 0.511, "Contact is at the wrong time (time=%)", getTime());
            contactHasOccured = true;
        }
    }
    
    void actionsAfterSolve() override
    {
        getHGrid()->info();
        particleHandler.write(std::cout);
        logger.assert_always(mathsFunc::isEqual(particleHandler.getObject(0)->getVelocity().X, -.94, 1e-2), "velocity"
                " not equal to restitution coefficient");
        logger.assert_always(mathsFunc::isEqual(particleHandler.getObject(1)->getVelocity().X, .94, 1e-2), "velocity"
                " not equal to restitution coefficient");
        logger.assert_always(contactHasOccured, "The contact has not occured");
    }

private:
    bool contactHasOccured;
};

int main(int argc, char* argv[])
{
    SphericalSuperQuadricCollision problem;
    problem.setName("SphericalSuperQuadricCollision");
    problem.setSaveCount(500);
    problem.solve();
    return 0;
}
