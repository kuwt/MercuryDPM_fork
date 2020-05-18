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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Boundaries/PeriodicBoundary.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/IntersectionOfWalls.h"
#include "Boundaries/SubcriticalMaserBoundary.h"
#include "Chute.h"
using namespace std;

class ChuteWithWedge : public Chute
{
public:
    ChuteWithWedge()
    {
        //Set and add the particle-species.
        LinearViscoelasticFrictionSpecies species;
        species.setDensity(6 / constants::pi);
        species.setCollisionTimeAndRestitutionCoefficient(5e-3, 0.88, 1);
        species.setSlidingDissipation(species.getDissipation());
        species.setSlidingStiffness(2. / 7. * species.getStiffness());
        species.setSlidingFrictionCoefficient(0.5);
        speciesHandler.copyAndAddObject(species);
    
        // Chute properties
        setChuteLength(20);
        setChuteWidth(10);
        setInflowHeight(10);
        inflowLength = 20;
        setInflowParticleRadius(0.5);
        setFixedParticleRadius(getInflowParticleRadius());
        setRoughBottomType(MULTILAYER);
        setChuteAngleAndMagnitudeOfGravity(27.0, 1);
    
        //make chute periodic in y-direction:
        makeChutePeriodic();
    
        //Set time-related properties
        setTimeMax(10);
        setTimeStep(1.0e-4); //(1/50th of the collision time)
        // Write 100 output files in total.
        setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(50, getTimeMax(), getTimeStep()));
    }
    
    void setupInitialConditions() override
    {
        Chute::setupInitialConditions();
        extendBottom();
        addFlowParticlesCompactly();
        addWedge();
        
        //add a maser boundary to generate particles:
        SubcriticalMaserBoundary b0;
        b0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMin() + inflowLength);
        boundaryHandler.copyAndAddObject(b0);
    }
    
    void addWedge()
    {
        IntersectionOfWalls w0;
        w0.setSpecies(speciesHandler.getObject(0));
        // wedge at an angle (90 degrees blocks the chute)
        Mdouble angle = 15 * constants::pi / 180;
        Vec3D normalIntoWall = Vec3D(-1, 0, 0);
        Vec3D point = Vec3D(70, 0, 0);
        w0.addObject(normalIntoWall, point);
        
        normalIntoWall = Vec3D(std::sin(angle), std::cos(angle), 0);
        point = Vec3D(60, 5, 0);
        w0.addObject(normalIntoWall, point);
        
        normalIntoWall = Vec3D(std::sin(angle), -std::cos(angle), 0);
        point = Vec3D(60, 5, 0);;
        w0.addObject(normalIntoWall, point);
        wallHandler.copyAndAddObject(w0);
    }
    
    void extendBottom()
    {
        const unsigned numberBottomParticles = particleHandler.getSize();
        const double initialLength = getChuteLength();
        setChuteLength(20 * getChuteLength());
        for (unsigned int i = 0; i < numberBottomParticles; ++i)
        {
            BaseParticle* particle = particleHandler.getObject(i);
            particle->setPosition(particle->getPosition() + Vec3D(initialLength, 0, 0));
            while (particle->getPosition().X < getChuteLength())
            {
                particleHandler.copyAndAddObject(particle);
                particle->setPosition(particle->getPosition() + Vec3D(initialLength, 0, 0));
            }
        }
    }
    
    
    ///creates flow particles in the whole chute
     SphericalParticle createFlowParticle() override
    {
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
        p0.setPosition(Vec3D(random.getRandomNumber(getXMin() + 2.0 * p0.getRadius(), inflowLength),
                             random.getRandomNumber(getYMin() + 2.0 * p0.getRadius(), getYMax()),
                             random.getRandomNumber(getZMin() + 2.0 * p0.getRadius(), getInflowHeight())));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        return p0;
    }
    Mdouble inflowLength;
};

int main(int argc, char* argv[])
{
    ChuteWithWedge problem;
    problem.setName("ChuteWithWedge");
    problem.setSaveCount(2000);
    
    problem.readArguments(argc, argv);
    problem.solve();
}
