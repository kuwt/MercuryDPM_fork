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

#include <iostream>
#include <iomanip> 
#include <vector>

#include "HGridCountCalls.h"
#include "Boundaries/PeriodicBoundary.h"
#include "HGridOptimiser.h"
#include "Species/LinearViscoelasticSpecies.h"

class HgridSpeedTest : public HGridCountCalls3D
{
public:
    HgridSpeedTest()
    {
        species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    }

    void loadParticlesFromDataFile(std::string& file)
    {
        if (readDataFile(file.c_str(), 7))
            std::cout << "inlezen gelukt" << std::endl;
        else
        {
            std::cout << "inlezen mislukt" << std::endl;
            exit(1);
        }
        
        boundaryHandler.clear();
        PeriodicBoundary b0;
        b0.set(Vec3D(1, 0, 0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        if (getSystemDimensions() >= 2)
        {
            b0.set(Vec3D(0, 1, 0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(b0);
        }
        if (getSystemDimensions() >= 3)
        {
            b0.set(Vec3D(0, 0, 1), getZMin(), getZMax());
            boundaryHandler.copyAndAddObject(b0);
        }
    }
    
    void createParticlesWithPowerLawDistribution(double omega, double alpha, double nu, int d, int N)
    {
        double VP;
        if (d == 1)
            VP = 2.0;
        if (d == 2)
            VP = constants::pi;
        if (d == 3)
            VP = constants::pi * 4.0 / 3.0;
        
        double L = std::pow(N * VP * DistInt(1, omega, d, omega, alpha) / nu, 1.0 / d);
        std::cout << "L=" << L << std::endl;
        setXMin(0);
        setYMin(0);
        setZMin(0);
        setXMax(L);
        setYMax(L);
        setZMax(L);
        
        PeriodicBoundary b0;
        b0.set(Vec3D(1, 0, 0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        if (d >= 2)
        {
            b0.set(Vec3D(0, 1, 0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(b0);
        }
        if (d >= 3)
        {
            b0.set(Vec3D(0, 0, 1), getZMin(), getZMax());
            boundaryHandler.copyAndAddObject(b0);
        }
        
        double V = 0;
        
        particleHandler.setStorageCapacity(2 * N);
        SphericalParticle p0;
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        
        //Use at least particles with maximum and minimum size
        p0.setRadius(1.0);
        if (d == 1)
            p0.setPosition(Vec3D(random.getRandomNumber(0, getXMax()), 0.0, 0.0));
        if (d == 2)
            p0.setPosition(Vec3D(random.getRandomNumber(0, getXMax()), random.getRandomNumber(0, getYMax()), 0.0));
        if (d == 3)
            p0.setPosition(Vec3D(random.getRandomNumber(0, getXMax()), random.getRandomNumber(0, getYMax()), random.getRandomNumber(0, getZMax())));
        particleHandler.copyAndAddObject(p0);
        V += VP * std::pow(p0.getRadius(), d);
        p0.setRadius(omega);
        if (d == 1)
            p0.setPosition(Vec3D(random.getRandomNumber(0, getXMax()), 0.0, 0.0));
        if (d == 2)
            p0.setPosition(Vec3D(random.getRandomNumber(0, getXMax()), random.getRandomNumber(0, getYMax()), 0.0));
        if (d == 3)
            p0.setPosition(Vec3D(random.getRandomNumber(0, getXMax()), random.getRandomNumber(0, getYMax()), random.getRandomNumber(0, getZMax())));
        particleHandler.copyAndAddObject(p0);
        V += VP * std::pow(p0.getRadius(), d);
        
        for (int i = 2; i < N; i++)
        {
            p0.setRadius(RandomRadius(omega, alpha));
            if (d == 1)
                p0.setPosition(Vec3D(random.getRandomNumber(0, getXMax()), 0.0, 0.0));
            if (d == 2)
                p0.setPosition(Vec3D(random.getRandomNumber(0, getXMax()), random.getRandomNumber(0, getYMax()), 0.0));
            if (d == 3)
                p0.setPosition(Vec3D(random.getRandomNumber(0, getXMax()), random.getRandomNumber(0, getYMax()), random.getRandomNumber(0, getZMax())));
            particleHandler.copyAndAddObject(p0);
            V += VP * std::pow(p0.getRadius(), d);
        }
        std::cout << "nu=" << V / std::pow(L, d) << " V=" << V << " L=" << L << std::endl;
    }
    
    double RandomRadius(double omega, double alpha)
    {
        double rand = random.getRandomNumber(0, 1);
        if (alpha == -1)
        {
            return std::pow(omega, rand);
        }
        else
        {
            return std::pow(rand * (std::pow(omega, 1.0 + alpha) - 1.0) + 1.0, 1.0 / (1.0 + alpha));
        }
    }

    double DistInt(double s, double e, int d, double omega, double alpha)
    {
        if (omega == 1)
        {
            return 1;
        }
        double teller;
        double noemer;
        if (alpha == -1)
        {
            noemer = std::log(omega);
        }
        else
        {
            noemer = (std::pow(omega, 1.0 + alpha) - 1.0) / (1.0 + alpha);
        }

        if (alpha + d == -1)
        {
            teller = std::log(e) - std::log(s);
        }
        else
        {
            teller = (std::pow(e, 1.0 + alpha + d) - std::pow(s, 1.0 + alpha + d)) / (1.0 + alpha + d);
        }
        return teller / noemer;
    }

    void callHGridRoutine()
    {
        for (std::vector<BaseParticle*>::iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
        {
            computeInternalForces(*it);
        }
    }
    
    void start()
    {
        perpareCalls();
        checkAndDuplicatePeriodicParticles();
        hGridActionsBeforeTimeLoop();
        hGridActionsBeforeTimeStep();
        callHGridRoutine();
        removeDuplicatePeriodicParticles();
        displayCalls(1);
    }
    
    const std::vector<double>& getHGridCellSizes()
    {
        return getHGrid()->getCellSizes();

    }

public:
    LinearViscoelasticSpecies* species;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    HgridSpeedTest problem;
    HGridOptimiser HGridOpt;
    
    problem.setHGridMaxLevels(8);
    {
        problem.random.randomise();
        problem.createParticlesWithPowerLawDistribution(20, -3, 0.4, 2, 1e4);
    }
    {
        //std::string filename="testPacking.data";
        //problem.loadParticlesFromDataFile(filename); 
    }
    
    std::cout << "Insutu testing:" << std::endl;
    
    problem.start();
    
    std::cout << "Analytical result:" << std::endl;
    HGridOpt.initialise(problem, 100, 0);
    std::vector<double> hGridCellSizes = problem.getHGridCellSizes();
    hGridCellSizes.insert(hGridCellSizes.begin(), 2.0 * problem.particleHandler.getSmallestParticle()->getRadius());
    HGridOpt.calculateWork(hGridCellSizes, problem.getHGridMethod(), 1);
    
}
