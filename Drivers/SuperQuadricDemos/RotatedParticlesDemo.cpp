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

#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "Species/HertzianViscoelasticMindlinSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Mercury3D.h"

class EllipticalSuperQuadricCollision : public Mercury3D
{
    void setupInitialConditions() override
    {
        
        LinearViscoelasticSlidingFrictionSpecies species;
        species.setCollisionTimeAndRestitutionCoefficient(5e-3, 0.8, 1);
        species.setSlidingStiffness(species.getStiffness());
        species.setSlidingDissipation(species.getDissipation());
        species.setDensity(constants::pi / 6);
        species.setSlidingFrictionCoefficient(0.5);
        speciesHandler.copyAndAddObject(species);
        setMax(8.0, 8.0, 20);
        setMin(1.5, 1.5, 1.5);
        setGravity({0, 0, -1});
        
        SuperQuadricParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setAxesAndExponents(2.0, 0.25, 0.25, 1.0, 1.0);
        p0.setInertia();
        
        unsigned int numberOfParticles = 50;
        unsigned int failed = 0;
        while (numberOfParticles > 0 && failed < 100)
        {
            Vec3D dummy;
            dummy.X = random.getRandomNumber(getXMin(), getXMax());
            dummy.Y = random.getRandomNumber(getYMin(), getYMax());
            dummy.Z = random.getRandomNumber(getZMin(), getZMax());
            p0.setPosition(dummy);
            
            dummy.X = random.getRandomNumber(-2, 2);
            dummy.Y = random.getRandomNumber(-2, 2);
            dummy.Z = random.getRandomNumber(-2, 2);
            p0.setVelocity(dummy);
            
            dummy.X = random.getRandomNumber(-1, 1);
            dummy.Y = random.getRandomNumber(-1, 1);
            dummy.Z = random.getRandomNumber(-1, 1);
            p0.setOrientationViaNormal(dummy);
            if (checkParticleForInteraction(p0))
            {
                particleHandler.copyAndAddObject(p0);
                numberOfParticles--;
            }
            else
            {
                failed++;
            }
        }
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0, 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0, 10.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(0.0, 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(10.0, 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        
        setTimeStep(1e-4);
        logger(INFO, "time step %", getTimeStep());
        setTimeMax(10);
        write(std::cout, true);
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
    problem.setName("EllipticalSuperQuadricCollision");
    problem.setSaveCount(500);
    problem.setSuperquadricParticlesWriteVTK(true);
    problem.setWallsWriteVTK(FileType::ONE_FILE);
    problem.solve();
    return 0;
}
