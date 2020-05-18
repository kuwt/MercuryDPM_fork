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

#include <iostream>
#include <iomanip> 
#include <chrono>
#include <Species/LinearViscoelasticSpecies.h>

#include "Mercury3D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Math/ExtendedMath.h"

class ScalingTestInitialConditionsRelax : public Mercury3D
{
public:
    
    void setupInitialConditions() override {
        SphericalParticle p0;
        p0.setRadius(particleRadius);
        
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        setXMax(std::pow(2.0 * N * 4.0 / 3.0 * constants::pi * std::pow(p0.getRadius(), 3.0), 1.0 / 3.0));
        setYMax(std::pow(2.0 * N * 4.0 / 3.0 * constants::pi * std::pow(p0.getRadius(), 3.0), 1.0 / 3.0));
        setZMax(std::pow(2.0 * N * 4.0 / 3.0 * constants::pi * std::pow(p0.getRadius(), 3.0), 1.0 / 3.0));
        
        for (int i = 0; i < N; i++)
        {
            double posX = random.getRandomNumber(getXMin(), getXMax());
            double posY = random.getRandomNumber(getXMin(), getXMax());
            double posZ = random.getRandomNumber(getXMin(), getXMax());
            
            double velX = random.getRandomNumber(-initialVelocity, initialVelocity);
            double velY = random.getRandomNumber(-initialVelocity, initialVelocity);
            double velZ = random.getRandomNumber(-initialVelocity, initialVelocity);
            
            p0.setPosition(Vec3D(posX, posY, posZ));
            p0.setVelocity(Vec3D(velX, velY, velZ));
            particleHandler.copyAndAddObject(p0);
        }
        
        PeriodicBoundary b0;
        b0.set(Vec3D(1, 0, 0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0, 1, 0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0, 0, 1), getZMin(), getZMax());
        boundaryHandler.copyAndAddObject(b0);
    }
    
    bool continueSolve() const override {
        Mdouble kineticEnergy = 0;
        
        for (std::vector<BaseParticle*>::const_iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
        {
            kineticEnergy += .5 * (*it)->getMass() * (*it)->getVelocity().getLengthSquared();
        }
        if (kineticEnergy < 0.5 * particleHandler.getLastObject()->getMass() * pow(targetVelocity, 2) * N)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    
    int N;
    double initialVelocity;
    double targetVelocity;
    double particleRadius;
};

class ScalingTestInitialConditionsEquilibrize : public Mercury3D
{
public:
    explicit ScalingTestInitialConditionsEquilibrize(Mercury3D& other)
            : DPMBase(other), Mercury3D(other)
    {
        
    }
};

class ScalingTestRun : public Mercury3D
{
public:
    explicit ScalingTestRun(Mercury3D& other)
            : DPMBase(other), Mercury3D(other)
    {
        multiplicationFactor_ = 1;
    }
    
    void setMultiplicationFactor(int multiplicationFactor)
    {
        multiplicationFactor_ = multiplicationFactor;
    }
    
    void setupInitialConditions() override {
        for (int j = 0; j < multiplicationFactor_ - 1; j++)
        {
            if (j % 3 == 0)
            {
                unsigned int N = particleHandler.getNumberOfObjects();
                for (int i = 0; i < N; i++)
                {
                    particleHandler.copyAndAddObject(particleHandler.getObject(i));
                    particleHandler.getLastObject()->move(Vec3D(getXMax(), 0.0, 0.0));
                }
                setXMax(2.0 * getXMax());
                dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(0))->moveRight(getXMax());
            }
            else if (j % 3 == 1)
            {
                
                unsigned int N = particleHandler.getNumberOfObjects();
                for (int i = 0; i < N; i++)
                {
                    particleHandler.copyAndAddObject(particleHandler.getObject(i));
                    particleHandler.getLastObject()->move(Vec3D(0.0, getYMax(), 0.0));
                }
                setYMax(2.0 * getYMax());
                dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(1))->moveRight(getYMax());
            }
            else if (j % 3 == 2)
            {
                
                unsigned int N = particleHandler.getNumberOfObjects();
                for (int i = 0; i < N; i++)
                {
                    particleHandler.copyAndAddObject(particleHandler.getObject(i));
                    particleHandler.getLastObject()->move(Vec3D(0.0, 0.0, getZMax()));
                }
                setZMax(2.0 * getZMax());
                dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(2))->moveRight(getZMax());
            }
        }
    }
    
    int multiplicationFactor_;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    ScalingTestInitialConditionsRelax scalingTestInitialConditionsRelax;
    auto species=scalingTestInitialConditionsRelax.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    scalingTestInitialConditionsRelax.N = 100;
    scalingTestInitialConditionsRelax.initialVelocity = 2.0;
    scalingTestInitialConditionsRelax.targetVelocity = 1.0;
    scalingTestInitialConditionsRelax.particleRadius = 0.5;
    scalingTestInitialConditionsRelax.setName("ScalingTestInitialConditionsRelax");
    species->setDensity(constants::pi / 6.0);
    scalingTestInitialConditionsRelax.setSaveCount(100);
    species->setStiffness(2e5);
    species->setDissipation(150);
    scalingTestInitialConditionsRelax.setTimeStep(1e-4);
    scalingTestInitialConditionsRelax.setTimeMax(1);
    scalingTestInitialConditionsRelax.setGravity(Vec3D(0.0, 0.0, 0.0));
    scalingTestInitialConditionsRelax.setHGridMaxLevels(1);
    scalingTestInitialConditionsRelax.solve();
    
    ScalingTestInitialConditionsEquilibrize scalingTestInitialConditionsEquilibrize(scalingTestInitialConditionsRelax);
    species=scalingTestInitialConditionsEquilibrize.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    scalingTestInitialConditionsEquilibrize.setName("ScalingTestInitialConditionsEquilibrize");
    species->setDissipation(0.0);
    scalingTestInitialConditionsEquilibrize.setSaveCount(100);
    scalingTestInitialConditionsEquilibrize.solve();
    
    for (unsigned int i = 1; i <= 10; i++)
    {
        ScalingTestRun scalingTestRun(scalingTestInitialConditionsEquilibrize);
        std::string name = "ScalingTestRun";
        name += std::to_string(i);
        scalingTestRun.setName(name);
        scalingTestRun.setMultiplicationFactor(i);
        scalingTestRun.setSaveCount(100);
        auto start = std::chrono::steady_clock::now();
        scalingTestRun.solve();
        auto end = std::chrono::steady_clock::now();
        auto diff = end - start;
        std::cout << "N="<<scalingTestRun.particleHandler.getNumberOfObjects()<<" T=" << std::scientific << std::setprecision(5) << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;
    }
    
}

