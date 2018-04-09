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

#include "Mercury2D.h"
#include "StatisticsVector.h"
#include "Boundaries/LeesEdwardsBoundary.h"
#include "Particles/BaseParticle.h"
#include "Species/LinearViscoelasticSpecies.h"

template <StatType T> class Relax : public StatisticsVector<T>, public Mercury2D
{
public:
    
    void setupInitialConditions()
    {
        Mdouble tc=1e-3;
        Mdouble r=0.8;
        Mdouble rho=2000;
        Mdouble velocity=0.1;
        Mdouble tmax=10;

        setName("Relax");
        setSystemDimensions(2);
        species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(rho);

        if(readDataFile("Startpos.dat",8))
                std::cout <<"inlezen gelukt"<< std::endl;
        else
        {
                std::cout <<"inlezen mislukt"<< std::endl;
                exit(1);
        }
        
        Mdouble minRadius=particleHandler.getSmallestParticle()->getRadius();
        Mdouble minMass=species->getMassFromRadius(minRadius);
        species->setStiffness(helpers::computeKFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass (tc, r, 0.5*minMass));
        species->setDissipation(helpers::computeDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass (tc, r, 0.5*minMass));

        setGravity(Vec3D(0.0,0.0,0.0));
        setTimeMax(tmax);
        
        setTimeStep(0.02*tc);
        getDataFile().setSaveCount(5000);
        getFStatFile().setFileType(FileType::NO_FILE);
        getEneFile().setSaveCount(500);
        getRestartFile().setSaveCount(5000);
        getStatFile().setSaveCount(500);
        
        setHGridMaxLevels(1);
        
        LR=boundaryHandler.copyAndAddObject(PeriodicBoundary());
        LR->set(Vec3D(1, 0, 0), getXMin(), getXMax());
        UD=boundaryHandler.copyAndAddObject(PeriodicBoundary());
        UD->set(Vec3D(0, 1, 0), getYMin(), getYMax());

    }

    
private:
    PeriodicBoundary* LR;
    PeriodicBoundary* UD;
    LeesEdwardsBoundary* leesEdwardsBoundary;
    LinearViscoelasticSpecies* species;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Relax<Z> relax;
    relax.setDoPeriodicWalls(false);
    relax.setDoTimeAverage(false);
    relax.solve();    
}
