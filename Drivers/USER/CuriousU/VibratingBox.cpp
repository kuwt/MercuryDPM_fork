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

#include <Mercury3D.h>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Boundaries/PeriodicBoundary.h>
using constants::pi;

/*
* This is our problem description. Everything is set up here.
* We inherit from Mercury3D, since this gives full flexibility.
* For more predefined problems (for instance, chutes), please see the
* documentation.
*/
class Drum : public Mercury3D
{
public:
    /* We define our own 'setupInitialConditions' function here,
     * which defines all the specifics about our simulation here.
     */
    void setupInitialConditions() override
    {
        //Particle radiusSmall
        Mdouble radiusSmall = 0.01;
        Mdouble radiusLarge = 0.02;
        //Vibration amplitude and frequency
        Mdouble amplitude = 1.0*radiusSmall;
        Mdouble frequency = 20; //Hz

        //The first step: set any properties which are always true for your system.
        // (for instance, if gravity remains constant, set it here)setName("Drum");
        setGravity(Vec3D(0,0,-9.8));
        setTimeMax(30);
        setTimeStep(0.0002);
        //visualised length
        setXMax(0.2);
        //visualised height
        setZMax(0.4);
        //visualised width
        setYMax(0.1);
        setXMin(0);
        setZMin(0);
        setYMin(-0.1);

        //Now, decide what Species you need for your system.
        LinearViscoelasticSlidingFrictionSpecies species;
        species.setDensity(10000);
        Mdouble massSmall = species.getDensity() * 4.0 / 3.0 * pi * radiusSmall * radiusSmall * radiusSmall;
        species.setCollisionTimeAndRestitutionCoefficient(25*getTimeStep(),0.8,massSmall);
        species.setSlidingStiffness(2./7.*species.getStiffness());
        species.setSlidingDissipation(2./7.*species.getDissipation());
        species.setSlidingFrictionCoefficient(0.5);
        speciesHandler.copyAndAddObject(species);

        //Add your walls below, and don't forget to set the species!

        //place a wall on the ground
        InfiniteWall ground;
        ground.setSpecies(speciesHandler.getObject(0));
        Vec3D position = Vec3D(0,0,0);
        Vec3D normal = Vec3D(0,0,-1);
        ground.set(normal,position);
        //Vibration function
        std::function<Vec3D(double)> vibration = [amplitude,frequency] (double time)
        {
            return Vec3D(0,0,amplitude * std::sin(time * 2 * pi * frequency));
        };
        ground.setPrescribedPosition(vibration);
        wallHandler.copyAndAddObject(ground);

        //place a left wall
        PeriodicBoundary periodicBoundary;
        periodicBoundary.set(Vec3D(1,0,0),0,getXMax());
        boundaryHandler.copyAndAddObject(periodicBoundary);

        //Place particles into the box
        SphericalParticle particle;
        particle.setSpecies(speciesHandler.getObject(0));
        //Insert a set of particles (3 large and 57 small):
        while (particleHandler.getNumberOfObjects()<60)
        {
            if (particleHandler.getNumberOfObjects()<3) {
                particle.setRadius(radiusLarge*random.getRandomNumber(0.9,1.1));
            } else {
                particle.setRadius(radiusSmall*random.getRandomNumber(0.9,1.1));
            }
            //Insert a particle:
            do
            {
                Mdouble x = random.getRandomNumber(0, getXMax());
                Mdouble z = random.getRandomNumber(0, getZMax());
                particle.setPosition(Vec3D(x,0,z));
            } while (!checkParticleForInteraction(particle));
            particleHandler.copyAndAddObject(particle);
        }
    }

};

int main(int argc, char **argv)
{
    Drum problem;
    problem.setName("VibratingBox");
    problem.setSaveCount(3*100);
    problem.wallHandler.setWriteVTK(FileType::ONE_FILE);
    problem.setParticlesWriteVTK(true);
    //problem.setParticlesWriteVTK(true);
    problem.setXBallsAdditionalArguments(" -v0 -solidf");
    problem.restartFile.setFileType(FileType::NO_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);

    //solve the system, the single particle will now bounce on the plate
    problem.solve(argc, argv);
    return 0;
}
