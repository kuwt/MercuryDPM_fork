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
#include "StatisticsVector.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include <sstream>
#include <iostream>
#include <cstdlib>
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"

template<StatType T> class pure_shear_untill_jamming : public StatisticsVector<T>, public Mercury2D
{
public:
    pure_shear_untill_jamming<T>()
            : StatisticsVector<T>(), Mercury2D()
    {
        species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    }

    void computeExternalForces(BaseParticle* CI)
    {
        /// Now add on gravity
        CI->addForce(getGravity() * CI->getMass());
        ///Finally walls
        //computeForcesDueToWalls(CI);
        ///Background friction
        CI->addForce(-CI->getVelocity() * 0.1 * species->getDissipation());
        CI->addTorque(-CI->getAngularVelocity() * 0.1 * species->getSlidingDissipation());
    }

    void setupInitialConditions()
    {
        if (readDataFile("Startpos.dat", 8))
            std::cout << "inlezen gelukt" << std::endl;
        else
        {
            std::cout << "inlezen mislukt" << std::endl;
            exit(1);
        }

        wallHandler.clear();
        boundaryHandler.clear();
        
        if (periodic)
        {
            PeriodicBoundary b0;
            b0.set(Vec3D(1, 0, 0), getXMin(), getXMax());
            boundaryHandler.copyAndAddObject(b0);
            b0.set(Vec3D(0, 1, 0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(b0);
        }
        else
        {
            InfiniteWall w0;
            w0.set(Vec3D(0.0, 1.0, 0.0), getYMax());
            wallHandler.copyAndAddObject(w0);
            w0.set(Vec3D(0.0, -1.0, 0.0), -getYMin());
            wallHandler.copyAndAddObject(w0);
            w0.set(Vec3D(1.0, 0.0, 0.0), getXMax());
            wallHandler.copyAndAddObject(w0);
            w0.set(Vec3D(-1.0, 0.0, 0.0), -getXMin());
            wallHandler.copyAndAddObject(w0);
        }
    }
    
    int cycles;
    bool periodic;

protected:
    void actionsAfterTimeStep()
    {
        double Xpressure = StatisticsVector<T>::getCGPoints()[0].NormalStress.XX + StatisticsVector<T>::getCGPoints()[0].TangentialStress.XX;
        double Ypressure = StatisticsVector<T>::getCGPoints()[0].NormalStress.YY + StatisticsVector<T>::getCGPoints()[0].TangentialStress.YY;

        if (Xpressure + Ypressure < 1e-3)
        {
            double fact = 1.0 + 1e-6;
            double xmax = getXMax() * fact;
            double ymax = getYMax() / fact;

            setXMax(xmax);
            setYMax(ymax);
            if (periodic)
            {
                dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(0))->moveRight(getXMax());
                dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(1))->moveRight(getYMax());
            }
            else
            {
                wallHandler.getObject(0)->setPosition(Vec3D(0.0, getYMax(), 0.0));
                wallHandler.getObject(2)->setPosition(Vec3D(getXMax(), 0.0, 0.0));
            }

        }
    }
public:
    LinearViscoelasticSlidingFrictionSpecies* species;
};

int main(int /*argc*/, char **/*argv[]*/)
{

    pure_shear_untill_jamming<Z> problem;

    problem.autoNumber();
    problem.setDoPeriodicWalls(false);
    problem.setDoTimeAverage(false);

    problem.species->setDensity(5);
    Mdouble mass = problem.species->getMassFromRadius(0.0064);
    problem.species->setStiffnessAndRestitutionCoefficient(1e6, 0.8, mass);
    //problem.setStiffnessAndRestitutionCoefficient(1e4,0.8,problem.getDensity()*constants::pi*pow(0.0037,2));
    problem.species->setSlidingFrictionCoefficient(0.5);

    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::ONE_FILE);
    problem.setGravity(Vec3D(0, 0, 0));

    problem.setTimeMax(10000.0);
    problem.setTimeStep(0.02 * helpers::computeCollisionTimeFromKAndDispAndEffectiveMass(problem.species->getStiffness(), problem.species->getDissipation(), 0.5 * mass));
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(10000, problem.getTimeMax(), problem.getTimeStep()));
    //problem.set_number_of_saves(10000);
    problem.getStatFile().setSaveCount(1);
    problem.setName("pure_shear_untill_jamming");
    problem.write(std::cout, false);
    problem.solve();
}

