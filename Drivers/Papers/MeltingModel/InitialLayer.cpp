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



#include "Mercury3D.h"
#include "Species/MeltableSpecies.h"
#include "Particles/MeltableParticle.h"
#include "Walls/InfiniteWall.h"
#include <random>
#include "Initialization.h"


class InitialBed : public Initialization {
    double timeMax = 1.0;
    Mdouble meanRadius = 55e-6/2.0; //27e-6; // in m
    Mdouble minRadius = 40e-6/2.0; //13.5e-6; //22.6e-6/2.0; //in m
    Mdouble maxRadius = 153e-6/2.0; //72e-6;  //153.4/2.0; in m
    Mdouble logMeanRadius = -10.4970; //-10.41;  //-10.3693;
    Mdouble logStdRadius = 0.6446; //0.2758;

public:
    void setupInitialConditions() override
    {
        setScaleFactor(1e4);

        setName("InitialLayer");
        removeOldFiles();
        setGravity(Vec3D(0.0, 0.0, -9.81));
        double xyScale = 60; // how many particles in x/y-direction
        double zScale = 4; // how many particles in z-direction
        setDomain(Vec3D(0,0,0),Vec3D(xyScale,xyScale,zScale)*2.0*meanRadius);

        auto species = speciesHandler.copyAndAddObject(MeltableSpecies());
        species->setDensity(density);
        species->setElasticModulus(elasticModulus);
        species->setPoissonRatio(poissonRatio);
        species->setDissipation(dissipation);
        species->setSolidHeatCapacity(solidHeatCapacity);
        species->setLiquidHeatCapacity(liquidHeatCapacity);

        setTimeStep(timeStep);
        setTimeMax(1.0);
        setSaveCount(timeMax/timeStep/200);
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::ONE_FILE);

        //Bottom Wall
        InfiniteWall w;
        w.setSpecies(species); //w0.setSpecies(speciesHandler.getObject(1));
        w.set(Vec3D(0.0, 0.0, -1.0), getMin());
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(-1.0, 0.0, 0.0), getMin());
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(1.0, 0.0, 0.0), getMax());
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(0.0, -1.0, 0.0), getMin());
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(0.0, 1.0, 0.0), getMax());
        wallHandler.copyAndAddObject(w);

        // insert particles
        std::mt19937 gen;
        //gen.seed(0);
        std::lognormal_distribution<> d(logMeanRadius, logStdRadius);
        //add particles until the volume to be added is zero
        MeltableParticle p;
        p.setSpecies(species);
        p.setRadius(maxRadius);
        p.setTemperature(meltingTemperature-temperatureInterval);
        Mdouble fillHeight = getZMin()+maxRadius;
        const double porosity = 0.55;
        double insertionVolume = porosity*getTotalVolume();
        while (insertionVolume > 0) {
            Mdouble x = random.getRandomNumber(getXMin()+p.getRadius(),getXMax()-p.getRadius());
            Mdouble y = random.getRandomNumber(getYMin()+p.getRadius(),getYMax()-p.getRadius());
            Mdouble z = random.getRandomNumber(getZMin()+p.getRadius(),fillHeight);
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                insertionVolume -= p.getVolume();
                do {
                    p.setRadius(d(gen));
                } while (p.getRadius() < minRadius || p.getRadius() > maxRadius);
                if (particleHandler.getNumberOfObjects() % 10 == 0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects() % 100 == 0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects() % 1000 == 0) std::cout << '\n';
            } else {
                fillHeight += 0.001*maxRadius;
            }
        }
        logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());

        write(std::cout,false);
    }

    bool continueSolve() const override {
        static int count = 0;
        if (count++==dataFile.getSaveCount()) {
            count=0;
            return getKineticEnergy() / getGravitationalEnergy() > 1e-5;
        }
        return true;
    }

    void printTime() const override {
        //logger(INFO, "time % timeMax % eneRatio % com %", getTime(), getTimeMax(), getKineticEnergy()/getGravitationalEnergy(),getCentreOfMass().Z);
        auto p = dynamic_cast<const MeltableParticle*>(particleHandler.getLargestParticle());
        logger.assert_debug(p,"MeltableParticle needed");
        logger(INFO, "time % timeMax % eneRatio % com % temp %", getTime(), getTimeMax(), getKineticEnergy()/getGravitationalEnergy(),getCentreOfMass().Z,p->getTemperature());
    }
};

int main(int argc UNUSED, char *argv[] UNUSED) {
    InitialBed problem;
    problem.solve();
    return 0;
}
