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
#include <initializer_list>
#include <vector>

#include "Mercury2D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Particles/BaseParticle.h"
#include "HGridOptimiser.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"

class HGridPaperFreeCooling : public Mercury2D
{
public:
    void printTime() const
    {
        std::cout << "t=" << std::setprecision(4) << std::left << std::setw(7) << getTime() << ", tmax=" << std::setprecision(4) << std::left << std::setw(7) << getTimeMax() << std::endl;
    }

    void setupInitialConditions()
    {
        int dimensions = 2;
        Mdouble density = 1000;
        Mdouble collisionTime = 1e-4;
        Mdouble restitutionCoefficient = 0.8;

        getFStatFile().setFileType(FileType::NO_FILE);
        getStatFile().setFileType(FileType::NO_FILE);

        //getDataFile().setFileType(FileType::NO_FILE);
        //getRestartFile().setFileType(FileType::NO_FILE);
        //getEneFile().setFileType(FileType::NO_FILE);

        getDataFile().setFileType(FileType::ONE_FILE);
        getRestartFile().setFileType(FileType::ONE_FILE);
        getEneFile().setFileType(FileType::ONE_FILE);
        getDataFile().setSaveCount(1000);
        getRestartFile().setSaveCount(1000);
        getEneFile().setSaveCount(100);

        setSystemDimensions(dimensions);
        setParticleDimensions(dimensions);

        setHGridDistribution(USER);
        setHGridMaxLevels(levelSizes.size());

        Mdouble effectiveMass;
        if (dimensions == 2)
            effectiveMass = 0.5 * density * constants::pi;
        else
            effectiveMass = 0.5 * density * constants::pi * 4.0 / 3.0;
        species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(density);
        species->setStiffness(helpers::computeKFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(collisionTime, restitutionCoefficient, effectiveMass));
        species->setDissipation(helpers::computeDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(collisionTime, restitutionCoefficient, effectiveMass));

        setTimeStep(0.02 * collisionTime);
        setTimeMax(0.1);
        setSaveCount(100);

        createParticlesWithPowerLawDistribution(dimensions);
        std::cout << "Finished creating particles" << std::endl;
    }
    
    void createParticlesWithPowerLawDistribution(int d)
    {
        double VP;
        switch (d)
        {
            case 1:
                VP = 2.0;
                break;
            case 2:
                VP = constants::pi;
                break;
            case 3:
                VP = constants::pi * 4.0 / 3.0;
                break;
            default:
                std::cerr << "Unknown dimension" << std::endl;
                exit(-1);
        }
        
        double L = std::pow(n * VP * DistInt(1, omega, d) / nu, 1.0 / d);
        setXMin(0);
        setYMin(0);
        setZMin(0);
        setXMax(L);
        setYMax(L);
        setZMax(L);
        
        if (periodic)
        {
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
        }
        else
        {
            InfiniteWall w0;
            w0.set(Vec3D(1.0, 0.0, 0.0), getXMax());
            wallHandler.copyAndAddObject(w0);
            w0.set(Vec3D(-1.0, 0.0, 0.0), -getXMin());
            wallHandler.copyAndAddObject(w0);
            if (d >= 2)
            {
                w0.set(Vec3D(0.0, 1.0, 0.0), getYMax());
                wallHandler.copyAndAddObject(w0);
                w0.set(Vec3D(0.0, -1.0, 0.0), -getYMin());
                wallHandler.copyAndAddObject(w0);
            }
            if (d >= 3)
            {
                w0.set(Vec3D(0.0, 0.0, 1.0), getZMax());
                wallHandler.copyAndAddObject(w0);
                w0.set(Vec3D(0.0, 0.0, -1.0), -getZMin());
                wallHandler.copyAndAddObject(w0);
            }
        }

        double V = 0;
        
        particleHandler.setStorageCapacity(2 * n);
        BaseParticle p0;
        p0.setSpecies(species);
        
        double absVel = 1000;

        for (int i = 0; i < n; i++)
        {
            p0.setRadius(RandomRadius(1.0 * (n - 1 - i) / (n - 1)));
            setRandVel(&p0, absVel, d);
            if (noOverlap)
            {
                do
                {
                    setRandPos(&p0, d);
                }
                while (!checkParticleForInteraction(p0));
            }
            else
            {
                setRandPos(&p0, d);
            }
            particleHandler.copyAndAddObject(p0);
            V += VP * std::pow(p0.getRadius(), d);
        }
    }
    
    void setRandPos(BaseParticle* P0, int d)
    {
        if (d == 1)
            P0->setPosition(Vec3D(random.getRandomNumber(0, getXMax()), 0.0, 0.0));
        if (d == 2)
            P0->setPosition(Vec3D(random.getRandomNumber(0, getXMax()), random.getRandomNumber(0, getYMax()), 0.0));
        if (d == 3)
            P0->setPosition(Vec3D(random.getRandomNumber(0, getXMax()), random.getRandomNumber(0, getYMax()), random.getRandomNumber(0, getZMax())));
    }

    void setRandVel(BaseParticle* P0, double absVel, int d)
    {
        if (d == 1)
            P0->setVelocity(Vec3D(random.getRandomNumber(-absVel, absVel), 0.0, 0.0));
        if (d == 2)
            P0->setVelocity(Vec3D(random.getRandomNumber(-absVel, absVel), random.getRandomNumber(-absVel, absVel), 0.0));
        if (d == 3)
            P0->setVelocity(Vec3D(random.getRandomNumber(-absVel, absVel), random.getRandomNumber(-absVel, absVel), random.getRandomNumber(-absVel, absVel)));
    }

    double RandomRadius(double rand = -1)
    {
        double rad;
        if (rand == -1)
        {
            rand = random.getRandomNumber(0, 1);
        }
        if (alpha == -1)
        {
            rad = std::pow(omega, rand);
        }
        else
        {
            rad = std::pow(rand * (std::pow(omega, 1.0 + alpha) - 1.0) + 1.0, 1.0 / (1.0 + alpha));
        }
        if (rad > omega)
            rad = omega;
        if (rad < 1.0)
            rad = 1.0;
        return rad;
    }
    
    double DistInt(double s, double e, int d)
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

    double userHGridCellSize(unsigned int level)
    {
        if (level >= levelSizes.size())
        {
            std::cerr << "In double HGridPaperFreeCooling::userHGridCellSize(unsigned int level) with level=" << level << std::endl;
        }
        return levelSizes[level];
    }

    unsigned int getHGridTargetNumberOfBuckets() const
    {
        return buckets;
    }

    double omega;
    int alpha;
    double nu;
    int n;
    bool periodic;
    bool noOverlap;
    double buckets;
    std::vector<double> levelSizes;
    LinearViscoelasticSpecies* species;
};

void createInitialCondition(double omega, int alpha, double nu, int N, bool periodic, bool noOverlap)
{
    HGridPaperFreeCooling problem;
    if (periodic && noOverlap)
    {
        std::cout << "Periodic walls and non overlapping insertion is not working (yet)" << std::endl;
        return;
    }
    if (noOverlap && nu > 0.6)
    {
        std::cout << "No overlap and large packing fraction not supported" << std::endl;
        return;
    }

    std::stringstream name;
    name << "RestartDataFreeCooling2DOmega" << omega << "Alpha" << alpha << "Nu" << nu << "N" << std::scientific << std::setprecision(0) << static_cast<double>(N);
    if (periodic)
        name << "Periodic";
    else
        name << "Walls";
    if (noOverlap)
        name << "NoOverlap";
    problem.setName(name.str());

    problem.omega = omega;
    problem.alpha = alpha;
    problem.nu = nu;
    problem.n = N;
    problem.periodic = periodic;
    problem.noOverlap = noOverlap;
    problem.buckets = 10 * N;

    //Optima for paramter K=0.3
    if (omega == 10 && alpha == -3 && nu == 0.4)
    {
        problem.levelSizes=
        {   4.0320, 7.8789, 14.4788};
    }
    else if( omega == 10 && alpha == -3 && nu == 0.6)
    {
        problem.levelSizes=
        {   3.7566, 6.8967, 12.2257};
    }
    else if( omega == 10 && alpha == -3 && nu == 0.8)
    {
        problem.levelSizes=
        {   3.5787, 6.2801, 10.7121, 17.2708};
    }
    else if (omega == 20 && alpha == -3 && nu == 0.4)
    {
        problem.levelSizes=
        {   4.3258, 9.1311, 18.6154, 35.4842};
    }
    else if(omega == 20 && alpha == -3 && nu == 0.6)
    {
        problem.levelSizes=
        {   4.0116, 7.8832, 15.0510, 27.2159};
    }
    else if(omega == 20 && alpha == -3 && nu == 0.8)
    {
        problem.levelSizes=
        {   3.8119, 7.1428, 13.1061, 23.3635};
    }
    else
    {
        std::cout<<"Unknown configuration"<<std::endl;
        return;
    }
    problem.levelSizes.push_back(nextafter(2 * omega, std::numeric_limits<double>::max()));

    problem.solve();
}

int main(int argc UNUSED, char *argv[] UNUSED)
{
    /*createInitialCondition(10, -3, 0.4, 1e2, true, false);
     createInitialCondition(10, -3, 0.4, 1e2, false, false);
     createInitialCondition(10, -3, 0.8, 1e2, true, false);
     createInitialCondition(10, -3, 0.8, 1e2, false, false);
     createInitialCondition(20, -3, 0.4, 1e2, true, false);
     createInitialCondition(20, -3, 0.4, 1e2, false, false);
     createInitialCondition(20, -3, 0.8, 1e2, true, false);
     createInitialCondition(20, -3, 0.8, 1e2, false, false);*/

    std::vector<double> omega = {10, 20};
    std::vector<double> nu = {0.4, 0.8};
    std::vector<double> N = {1e2, 1e3, 1e4};
    for (auto itOmega = omega.begin(); itOmega != omega.end(); ++itOmega)
    {
        for (auto itNu = nu.begin(); itNu != nu.end(); ++itNu)
        {
            for (auto itN = N.begin(); itN != N.end(); ++itN)
            {
                createInitialCondition(*itOmega, -3, *itNu, *itN, true, true);
                createInitialCondition(*itOmega, -3, *itNu, *itN, false, true);
                createInitialCondition(*itOmega, -3, *itNu, *itN, true, false);
                createInitialCondition(*itOmega, -3, *itNu, *itN, false, false);
            }
        }
    }

    createInitialCondition(20, -3, 0.6, 1e3, true, true);
}
