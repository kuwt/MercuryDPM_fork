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

#include <Boundaries/CubeInsertionBoundary.h>
#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <MercuryTime.h>
#include <Walls/InfiniteWall.h>

class DPM : public Mercury3D {
    void printTime() const override {
        logger(INFO,"t % ene % pressure %",getTime(),getKineticEnergy()/getElasticEnergy(),getStaticStress().trace()/3.0);
    }

    bool continueSolve() const override {
        return getNumberOfTimeSteps()<1000 || getKineticEnergy()/getElasticEnergy()>1e-5;
    }

public:
    void computeAllForces() override {
        Mercury3D::computeAllForces();
    }

    void computeAllForcesNoHGrid(){
        //Resetting all forces on both particles and walls to zero
        for (BaseParticle* const p : particleHandler)
        {
            p->setForce(Vec3D(0.0, 0.0, 0.0));
            p->setTorque(Vec3D(0.0, 0.0, 0.0));
        }
        for (BaseWall* const w : wallHandler)
        {
            w->setForce(Vec3D(0.0, 0.0, 0.0));
            w->setTorque(Vec3D(0.0, 0.0, 0.0));
        }
        for (auto i : interactionHandler) {
            auto w = dynamic_cast<BaseWall*>(i->getI());
            if (w == nullptr) {
                DPMBase::computeInternalForce((BaseParticle*)i->getP(), (BaseParticle*)i->getI());
            } else {
                DPMBase::computeForcesDueToWalls((BaseParticle*)i->getP(), w);
            }
        }
        for (auto p : particleHandler) {
            computeExternalForces(p);
        }
    }

//    void computeInternalForce(BaseParticle* const p, BaseParticle* const q){
//        DPMBase::computeInternalForce(p, q);
//    }
//
//    void computeForcesDueToWalls(BaseParticle* p, BaseWall* w) {
//        DPMBase::computeForcesDueToWalls(p, w)
//    }
};

int main()
{
    // create a simulation
    DPM dpm;
    // restart from initial condition; create them if they dont exist.
    if (!dpm.readRestartFile("InteractionHandlerSpeedTestInit")) {
        logger(ERROR,"Creating restart file first");
        // create a initial simulation
        unsigned n = 100000;
        DPM dpm0;
        dpm0.setName("InteractionHandlerSpeedTestInit");
        dpm0.setDomain(Vec3D(0, 0, 0), cbrt(n) * 0.95 * Vec3D(1, 1, 1));
        //dpm0.setGravity(Vec3D(0,0,-1));
        dpm0.setTimeStep(1e-4);
        dpm0.setTimeMax(10);
        dpm0.setSaveCount(50);
        // create species
        auto species = dpm0.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(6.0/constants::pi);
        species->setCollisionTimeAndRestitutionCoefficient(10 * dpm0.getTimeStep(), 0.1, 1.0);
        // create walls
        InfiniteWall wall(species);
        wall.set(Vec3D(0,0,1), dpm0.getMax());
        dpm0.wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(0,0,-1), dpm0.getMin());
        dpm0.wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(0,1,0), dpm0.getMax());
        dpm0.wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(0,-1,0), dpm0.getMin());
        dpm0.wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(1,0,0), dpm0.getMax());
        dpm0.wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(-1,0,0), dpm0.getMin());
        dpm0.wallHandler.copyAndAddObject(wall);
        // create particles
        SphericalParticle particle(species);
        CubeInsertionBoundary insertionBoundary;
        insertionBoundary.set(particle, 1e6, dpm0.getMin(), dpm0.getMax(), {0, 0, 0}, {0, 0, 0}, .47, .53);
        insertionBoundary.setInitialVolume((n-0.5)*constants::pi/6.);
        insertionBoundary.checkBoundaryBeforeTimeStep(&dpm0);
        logger(INFO, "Number of particles: %", dpm0.particleHandler.getSize());
        //solve
        dpm0.solve();
        if (!dpm.readRestartFile("InteractionHandlerSpeedTestInit")) {
            logger(ERROR, "Restart file could not be created");
        }
    }

    //rename simulation
    dpm.setName("InteractionHandlerSpeedTest");
    dpm.setTime(0);
    //simulate 0 time steps
    dpm.setTimeMax(0);
    dpm.setHGridMaxLevels(1);
    dpm.setFileType(FileType::NO_FILE);
    dpm.solve();
    dpm.interactionHandler.clear();
    dpm.computeAllForces();
    Time timer;
    int repetitions = 10;
    logger(INFO,"% interactions",dpm.interactionHandler.getNumberOfObjects());
    dpm.hGridInfo(std::cout);


    //timings
    timer.tic();
    for (int i = 0; i < repetitions; ++i) {
        dpm.interactionHandler.clear();
        dpm.computeAllForces();
    }
    logger(INFO, "Time to simulate contact detection and force computation for % non-existing contacts: % s",dpm.interactionHandler.getNumberOfObjects(), timer.toctic());

    for (int i = 0; i < repetitions; ++i) {
        dpm.computeAllForces();
    }
    logger(INFO, "Time to simulate contact detection and force computation for % existing contacts: % s",dpm.interactionHandler.getNumberOfObjects(), timer.toctic());

    for (int i = 0; i < repetitions; ++i) {
        dpm.computeAllForcesNoHGrid();
    }
    logger(INFO, "Time to simulate force computation for % existing contacts: % s",dpm.interactionHandler.getNumberOfObjects(), timer.toctic());

    return 0;
}

