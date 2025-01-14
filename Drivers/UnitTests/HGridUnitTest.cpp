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

#include <Species/LinearViscoelasticSpecies.h>
#include "Mercury3D.h"
#include "DPMBase.h"
#include "Boundaries/PeriodicBoundary.h"
class MD_demo : public DPMBase
{
public:

    MD_demo()
    {
        species = new LinearViscoelasticSpecies;
        speciesHandler.addObject(species);
    }

    void setupInitialConditions() override {
        if (particleHandler.getNumberOfObjects() != N)
        {
            particleHandler.clear();
            boundaryHandler.clear();
            
            const double VP = constants::pi * 4.0 / 3.0;
            L = pow(N * VP * DistInt(1, omega) / nu, 1.0 / 3.0);
            
            setXMin(0);
            setYMin(0);
            setZMin(0);
            setXMax(L);
            setYMax(L);
            setZMax(L);
            
            PeriodicBoundary b0;
            b0.set(Vec3D(1, 0, 0), getXMin(), getXMax());
            boundaryHandler.copyAndAddObject(b0);
            b0.set(Vec3D(0, 1, 0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(b0);
            b0.set(Vec3D(0, 0, 1), getZMin(), getZMax());
            boundaryHandler.copyAndAddObject(b0);
            
            particleHandler.setStorageCapacity(2 * N);
            SphericalParticle p0;
            p0.setSpecies(speciesHandler.getObject(0));
            p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
            
            double V = 0;
            
            //Use at least particles with maximum and minimum size
            p0.setRadius(1.0);
            Vec3D position;
            position.Z = random.getRandomNumber(0, getZMax());
            position.Y = random.getRandomNumber(0, getYMax());
            position.X = random.getRandomNumber(0, getXMax());
            p0.setPosition(position);
            particleHandler.copyAndAddObject(p0);
            V += VP * pow(p0.getRadius(), 3);
            
            p0.setRadius(omega);
            position.Z = random.getRandomNumber(0, getZMax());
            position.Y = random.getRandomNumber(0, getYMax());
            position.X = random.getRandomNumber(0, getXMax());
            p0.setPosition(position);
            particleHandler.copyAndAddObject(p0);
            V += VP * pow(p0.getRadius(), 3);

            //For other particles use a random distribution
            for (unsigned int i = 2; i < N; i++)
            {
                p0.setRadius(RandomRadius());
                position.Z = random.getRandomNumber(0, getZMax());
                position.Y = random.getRandomNumber(0, getYMax());
                position.X = random.getRandomNumber(0, getXMax());
                p0.setPosition(position);
                particleHandler.copyAndAddObject(p0);
                V += VP * pow(p0.getRadius(), 3);
            }
        }
    }
    
    double RandomRadius()
    {
        double rand = random.getRandomNumber(0, 1);
        if (alpha == -1)
        {
            return pow(omega, rand);
        }
        else
        {
            return pow(rand * (pow(omega, 1.0 + alpha) - 1.0) + 1.0, 1.0 / (1.0 + alpha));
        }
    }
    
    double DistInt(const double s, const double e)
    {
        if (omega == 1)
        {
            return 1;
        }
        double numerator;
        double denominator;
        if (alpha == -1)
        {
            denominator = log(omega);
        }
        else
        {
            denominator = (pow(omega, 1.0 + alpha) - 1.0) / (1.0 + alpha);
        }
        
        if (alpha == -4)
        {
            numerator = log(e) - log(s);
        }
        else
        {
            numerator = (pow(e, 4.0 + alpha) - pow(s, 4.0 + alpha)) / (4.0 + alpha);
        }
        return numerator / denominator;
    }
    
    double L;
    double omega;
    double alpha;
    double nu;
    unsigned int N;
    LinearViscoelasticSpecies* species;
};

class HGrid_demo : public Mercury3D
{
public:
    explicit HGrid_demo(DPMBase& other)
            : DPMBase(other)
    {
    }
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    MD_demo MD_problem;
    MD_problem.species->setDensity(2000);

    MD_problem.omega = 40;
    MD_problem.alpha = -2;
    MD_problem.nu = 0.7;
    MD_problem.N = 1000;
    MD_problem.setSystemDimensions(3);
    MD_problem.setParticleDimensions(3);
    MD_problem.setName("HGridUnitTest_MD");
    MD_problem.setTimeStep(0.1);
    MD_problem.setTimeMax(0.09);
    MD_problem.setSaveCount(1);
    MD_problem.setupInitialConditions();

    HGrid_demo HGrid_problem1(MD_problem);
    HGrid_problem1.setHGridMethod(TOPDOWN);
    HGrid_problem1.setHGridMaxLevels(3);
    HGrid_problem1.setHGridDistribution(EXPONENTIAL);
    HGrid_problem1.setName("HGridUnitTest_HGrid1");
    
    HGrid_demo HGrid_problem2(MD_problem);
    HGrid_problem2.setHGridMethod(BOTTOMUP);
    HGrid_problem2.setHGridMaxLevels(8);
    HGrid_problem2.setHGridDistribution(LINEAR);
    HGrid_problem2.setName("HGridUnitTest_HGrid2");

    logger(INFO, "Solving the MD problem");
    MD_problem.solve();
    logger(INFO, "Solving the first HGrid problem");
    HGrid_problem1.solve();
    logger(INFO, "Solving the second HGrid problem");
    HGrid_problem2.solve();
    
    // Check the particles are in the same place for all three problem i.e. no HGrid and two different HGrid settings
    const double tolerance = 1e-10;
    std::vector<BaseParticle*>::iterator hGrid1It = HGrid_problem1.particleHandler.begin();
    std::vector<BaseParticle*>::iterator hGrid2It = HGrid_problem2.particleHandler.begin();
    for(BaseParticle* particle : MD_problem.particleHandler)
    {
        if (!(particle->getPosition().isEqualTo((*hGrid1It)->getPosition(), tolerance)))
        {
            logger(ERROR, "position of particle in hGrid 1 is not equal to the position without hGrid");
        }
        if (!(particle->getPosition().isEqualTo((*hGrid2It)->getPosition(), tolerance)))
        {
            logger(ERROR, "position of particle in hGrid 2 is not equal to the position without hGrid");
        }
        ++hGrid1It;
        ++hGrid2It;
    }
}
