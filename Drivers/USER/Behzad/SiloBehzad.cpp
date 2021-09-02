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

#include "Diameter.h"
#include <Mercury3D.h>
#include <Species/HertzianViscoelasticMindlinRollingTorsionSpecies.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/IntersectionOfWalls.h>
#include <Boundaries/InsertionBoundary.h>
#include <Boundaries/CubeDeletionBoundary.h>
#include <cstring>

//scale up factor (set to one for real simulations, set higher for running a demonstration)
Mdouble scaleup = 1;
const Mdouble fillTime_ = 1./scaleup;
const Mdouble settlingTime = 0.2;

/*
 * Here, we create a special insertion boundary that is similar to the cubeInsertionBoundary:
 * - a limit to the amount of particles is imposed by adding particles only if the amount of particles in the system is below NFinal/fillTime
 * - at most one particle is inserted per time step
 * - the radius of the inserted particle is taken from the diameter array; the i-th particle is of radius 0.5*diameter[i].
 */
class BehzadInsertionBoundary : public InsertionBoundary
{
    BaseParticle* generateParticle(RNG& random)
    {
        DPMBase& dpm = *getHandler()->getDPMBase();
        BaseParticle* P = getParticleToCopy()->copy();

        P->setRadius(0.5 * diameter[dpm.particleHandler.getSize()%diameter.size()]*scaleup);

        return P;
    }

    void placeParticle(BaseParticle* p, RNG& random)
    {
        DPMBase& dpm = *getHandler()->getDPMBase();
        Vec3D min(dpm.getXMin()+p->getRadius(),dpm.getYMin(),dpm.getZMax());
        Vec3D max(dpm.getXMax()-p->getRadius(),dpm.getYMax(),1.3*dpm.getZMax());
        Vec3D pos;
        pos.X = random.getRandomNumber(min.X, max.X);
        pos.Y = random.getRandomNumber(min.Y, max.Y);
        pos.Z = random.getRandomNumber(min.Z, max.Z);
        p->setPosition(pos);
    }

    void checkBoundaryBeforeTimeStep(DPMBase* md) {
        unsigned int failed = 0;
        DPMBase& dpm = *getHandler()->getDPMBase();
        unsigned long size = std::min(1.0,dpm.getTime()/fillTime_)*33000/scaleup/scaleup/scaleup;
        if (dpm.getTime()==0) logger(INFO,"Number of particles to insert %", 33000/scaleup/scaleup/scaleup);
        //setting a max rate of 1/40 per time step
        //try max_failed times to find new insertable particle
        while (failed <= maxFailed_ && dpm.particleHandler.getNumberOfObjects()<size) {
            BaseParticle* p0 = generateParticle(md->random);
            p0->setHandler(&md->particleHandler);

            if (md->checkParticleForInteraction(*p0)) {
                BaseParticle *added = md->particleHandler.copyAndAddObject(p0);
                added->setSpecies(md->speciesHandler.getObject(
                        added->getIndSpecies()
                ));
                failed = 0;
                ++numberOfParticlesInserted_;
            } else {
                failed++;
            }
        }
    }

public:
    BehzadInsertionBoundary* copy() const {
        return new BehzadInsertionBoundary(*this);
    }

    std::string getName() const {
        return "BehzadInsertionBoundary";
    }
};

/*
 * This is the definition of the dpm problem:
 * - in setupInitialConditions, the initial walls (side, base) and boundary conditions (periodicity, insertion) are set;
 * - in actionsBeforeTimeStep, two changes to the setup are made, (1) after filling and (2) after settling of the particles.
 */
class Silo final : public Mercury3D
{
    //this radius is used to set the dimensions of the silo (height, width, length) and orifice (width), as well as the restitution coefficient
    const Mdouble radius_ = 1.5e-3;

public:

    Mdouble muRolling = 0;

    void setupInitialConditions() override {

        //set general properties
        setMin(Vec3D(-50, -5, 0) * radius_);
        setMax(Vec3D(50, 5, 200) * radius_);
        setGravity(Vec3D(0, 0, -9.8));

        //define material/contact properties
        HertzianViscoelasticMindlinRollingTorsionSpecies s;
        s.setDensity(950);

        Mdouble poisson = 0.4;
        Mdouble shearModulus = 1e8;
        Mdouble restitutionCoeff = 1e-4;
        s.setEffectiveElasticModulusAndRestitutionCoefficient(
                shearModulus * (3 * poisson + 2 * shearModulus) / (poisson + shearModulus), restitutionCoeff);

        s.setSlidingFrictionCoefficient(.577);
        s.setRollingFrictionCoefficient(muRolling);
        s.setTorsionFrictionCoefficient(0);
        s.setEffectiveShearModulus(shearModulus);
        ParticleSpecies* particleSpecies = speciesHandler.copyAndAddObject(s);
        ParticleSpecies* wallSpecies = speciesHandler.copyAndAddObject(s);
        dynamic_cast<HertzianViscoelasticMindlinRollingTorsionMixedSpecies*>(speciesHandler.getMixedObject(0,1))
                ->setSlidingFrictionCoefficient(.5);

        //periodic in y-direction
        PeriodicBoundary b;
        b.set(Vec3D(0, 1, 0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b);

        //walls in x-direction
        wallHandler.copyAndAddObject(InfiniteWall(Vec3D(+1, 0, 0), getMax(), wallSpecies));
        wallHandler.copyAndAddObject(InfiniteWall(Vec3D(-1, 0, 0), getMin(), wallSpecies));

        SphericalParticle p;
        p.setSpecies(particleSpecies);

        BehzadInsertionBoundary c; //delete is done in boundaryHandler
        c.set(&p,0);
        boundaryHandler.copyAndAddObject(c);

        //add walls with opening in z-direction
        IntersectionOfWalls wLeft;
        wLeft.setSpecies(wallSpecies);
        wLeft.addObject(Vec3D(0, 0, -1), getMin());
        wLeft.addObject(Vec3D(-1, 0, 0), Vec3D(-10, 0, 0) * radius_);
        wallHandler.copyAndAddObject(wLeft);
        IntersectionOfWalls wRight;
        wRight.setSpecies(wallSpecies);
        wRight.addObject(Vec3D(0, 0, -1), getMin());
        wRight.addObject(Vec3D(1, 0, 0), Vec3D(10, 0, 0) * radius_);
        wallHandler.copyAndAddObject(wRight);
        IntersectionOfWalls closure;
        closure.setSpecies(wallSpecies);
        closure.addObject(Vec3D(0, 0, -1), getMin());
        closure.addObject(Vec3D(-1, 0, 0), Vec3D(10, 0, 0) * radius_);
        closure.addObject(Vec3D(1, 0, 0), Vec3D(-10, 0, 0) * radius_);
        wallHandler.copyAndAddObject(closure);

        //time-stepping properties
        setTimeStep(1.4e-6);
        setTimeMax(10.0);

        //output properties
        setSaveCount(2143);
        restartFile.setSaveCount(21430);
        restartFile.setFileType(FileType::MULTIPLE_FILES);
        setWallsWriteVTK(FileType::ONE_FILE);
    }

    /*
     * This function modifies the output to the console. Currently, time, particle number, and the relaxation state is written out.
     */
    void printTime() const override
    {
        logger(INFO,"t %\tN %\tr %\tE %",getTime(),particleHandler.getNumberOfObjects(),round(helpers::getRealTime()),getKineticEnergy()/getElasticEnergy());
    }

    /*
     * This function makes two interventions in the code:
     * - when fillingTime is reached, we remove the insertionBoundary, remove all particles above zMax, and set the openTime to currentTime+settlingTime
     * - when openTime is reached, we remove the wall closing the orifice and add a deletion boundary that removes particles below -0.1*zMax
     * The user is informed of these events via logger messages to the console
     */
    void actionsBeforeTimeStep() override {
        static bool closed = true;
        static bool filling = true;
        static Mdouble openTime = 0;
        if (filling) {
            if (getTime()>=fillTime_) {
                logger(INFO, "Stop filling (N=%), remove extra particles",particleHandler.getSize());
                filling = false;
                openTime = getTime() + settlingTime;
                boundaryHandler.removeLastObject();
                for (int i = particleHandler.getSize() - 1; i != -1; --i) {
                    if (particleHandler.getObject(i)->getPosition().Z > getZMax()) {
                        particleHandler.removeObject(i);
                    }
                }
            }
        } else if (closed && getTime()>=openTime) {
            logger(INFO,"Opening outflow (N=%)",particleHandler.getSize());
            closed = false;
            //remove closure wall, allowing particles to flow out
            wallHandler.removeLastObject();
            //add deletion boundary, removing all particles below a height of -0.1*getZMax()
            CubeDeletionBoundary db;
            db.set(getMin()-Vec3D(0,0,getZMax()),getMax()-Vec3D(0,0,1.1*getZMax()));
            boundaryHandler.copyAndAddObject(db);
            writeRestartFile();
        }
    }

    // write new read-in options
    bool readNextArgument(int& i, int argc, char *argv[]) override
    {
        if (!strcmp(argv[i], "-murolling"))
        {
            muRolling = atof(argv[i + 1]);
            logger(INFO,"set muRolling=%",muRolling);
            setName(getName()+"MuR"+argv[i + 1]);
            logger(INFO,"set name=%",getName());
        }
        else {
            return Mercury3D::readNextArgument(i, argc, argv);
        } //if argv[i] is not found, check for commands in the BaseClass

        return true; //returns true if argv[i] is found
    }

};

int main(int argc, char** argv)
{
    Silo silo;
    silo.setName("SiloBehzad");
    silo.solve(argc, argv);
    return 0;
}
