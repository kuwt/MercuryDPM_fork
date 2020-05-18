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

#include "Mercury3D.h"
#include "Particles/SuperQuadricParticle.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Math/ExtendedMath.h"

class EllipticalSuperQuadricCollision : public Mercury3D
{
    void setupInitialConditions() override
    {
        
        LinearViscoelasticSpecies species;
        species.setStiffness(2e5);
        species.setDissipation(25);
        speciesHandler.copyAndAddObject(species);
        
        SuperQuadricParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setAxesAndExponents(2.0,1.0,1.0,1.0,1.0);
        p0.setInertia();
        SuperQuadricParticle p1 = *(p0.copy());
        
        p0.setPosition(Vec3D(0, 0.0, 0.0));
        p1.setPosition(Vec3D(2.5, 2.5, 0.0));
        p0.setVelocity(Vec3D(1, 1.0, 0.0));
        p1.setVelocity(Vec3D(-1, -1.0, 0.0));
        p0.setOrientationViaNormal({0, 0, 1});
        p1.setOrientationViaNormal({0,0,1});
        logger(INFO, "Euler angles p1: %",p1.getOrientation().getEuler());
        particleHandler.copyAndAddObject(p0);
        particleHandler.copyAndAddObject(p1);
        logger(INFO, "interaction radius p0: %", p0.getMaxInteractionRadius());
        
        
        setTimeStep(species.getCollisionTime(1) / 50);
        setTimeMax(0.6);
        setMin(-10, -1, -1);
        setMax(10, 1, 1);
    }
    
    //check contact time
    void actionsAfterTimeStep() override
    {
        
    }
    
    void actionsAfterSolve() override
    {
        getHGrid()->info();
        particleHandler.write(std::cout);
    }

private:
    
};

int main(int argc, char* argv[])
{
    EllipticalSuperQuadricCollision problem;
    problem.setName("EllipsoidSelfTest");
    problem.setSaveCount(50);
    // comment the  next to lines to turn on file output
    problem.setSuperquadricParticlesWriteVTK(true);
    problem.setFileType(FileType::NO_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.solve();
    return 0;
}
