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

#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "Mercury3D.h"
#include "Walls/Screw.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/InfiniteWall.h"

/**
 * This class defines all properties ot the Mixer in setupinitialcondititons, the the user can change then after class instantiation
 */
class NautaMixer : public Mercury3D
{
private:
    //The z-value of the tip of the cone that forms the outer mixer wall
    Mdouble coneTipHeight = 0;
    //The z-value of the tip of the cone that the screw rotates around
    Mdouble screwTipHeight = 0.5;
    //The z-value of the intersection of the cone and the base.
    Mdouble coneBaseHeight = 0.6;
    //The z-value of the tip of the cone that forms the base
    Mdouble baseTipHeight = 0.7;
    //The z-value of the top plate
    Mdouble topHeight = 2.5;

    //The distance from the z-axis of the intersection of the screw rotation axis and the top plate
    Mdouble screwTopRadius = 0.55;
    //The radius of the intersection of the cone and the top plate.
    Mdouble coneTopRadius = 0.7;

    //The radius of the screw
    Mdouble screwRadius = 0.12;
    //The radius of the screw
    Mdouble screwPitchNum = 9;
    //The thickness of the screw blade
    Mdouble screwThickness = 0.01;
    //The radius of the screw
    Mdouble screwCoreRadius = 0.025;

    //The rotation speed of the screw around the z-axis in rev/s
    Mdouble rotationRate = 0.2;
    //The (additional) rotation speed of the screw in rev/s of the screw
    Mdouble screwRotationRate = 6.0*rotationRate;

    //min particle radius
    Mdouble particleRadius = 0.025;
    //Difference between smallest and largest particle radius (uniform size distribution)
    Mdouble particlePolydispersity = 1.4;
    //maximal number of particles introduced
    unsigned maxParticleNumber = 20000;

public:
    void setupInitialConditions() override
    {
        removeOldFiles();
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::MULTIPLE_FILES);
        setDomain(Vec3D(-coneTopRadius, -coneTopRadius, coneTipHeight),
                  Vec3D(coneTopRadius, coneTopRadius, topHeight));
        setGravity({0,0,-9.8});

        addSpeciesAndSetTimeStepAndSaveCount();
        //setTimeMax(1 / rotationRate);

        addBaseWall();
        addConeWall();
        //addTopWall();
        addScrew();

        addParticles();
    }

//    void actionsAfterSolve() override {
//        //addParticlesAtWall();
//        addParticles();
//        forceWriteOutputFiles();
//    }

    //create a species for all walls and particles
    void addSpeciesAndSetTimeStepAndSaveCount() {
        LinearViscoelasticSlidingFrictionSpecies s;
        s.setHandler(&speciesHandler);
        s.setDensity(2000);
        Mdouble mass = s.getMassFromRadius(particleRadius);
        s.setCollisionTimeAndRestitutionCoefficient(0.003, 0.01, mass);
        s.setSlidingFrictionCoefficient(0.5);
        s.setSlidingStiffness(2.0/7.0*s.getStiffness());
        s.setSlidingDissipation(2.0/7.0*s.getDissipation());
        speciesHandler.copyAndAddObject(s);
        setTimeStep(0.1 * s.getCollisionTime(mass));
        logger(INFO,"Timestep %",getTimeStep());
        setSaveCount(10.0*s.getCollisionTime(mass)/getTimeStep());
        logger(INFO,"Saving every % s",getTimeStep()*dataFile.getSaveCount());
    }

    //define base wall and add to wallHandler
    void addBaseWall() {
        Mdouble coneBaseRadius = coneTopRadius * coneBaseHeight/topHeight;

        //two points on the base
        Vec3D tip = {0,0,baseTipHeight};
        Vec3D top = {coneBaseRadius,0,coneBaseHeight};
        //normal
        Vec3D normal = {top.Z-tip.Z,0,-top.X};

        AxisymmetricIntersectionOfWalls w;
        w.setAxis({0,0,1});
        w.addObject(normal,tip);
        w.setSpecies(speciesHandler.getLastObject());
        wallHandler.copyAndAddObject(w);
    }

    //define cone wall and add to wallHandler
    void addConeWall() {
        //two points on the cone
        Vec3D tip = {0,0,coneTipHeight};
        Vec3D top = {coneTopRadius,0,topHeight};
        //normal
        Vec3D normal = {top.Z-tip.Z,0,-top.X};

        AxisymmetricIntersectionOfWalls w;
        w.setAxis({0,0,1});
        w.addObject(normal,tip);
        w.setSpecies(speciesHandler.getLastObject());
        wallHandler.copyAndAddObject(w);
    }

    //define top wall and add to wallHandler
    void addTopWall() {
        InfiniteWall w;
        w.set({0,0,1},{0,0,topHeight});
        w.setSpecies(speciesHandler.getLastObject());
        wallHandler.copyAndAddObject(w);
    }

    //define screw  and add to wallHandler
    void addScrew()
    {
        //two points on the screw axis
        Vec3D tip = {0,0,screwTipHeight};
        Vec3D top = {screwTopRadius,0,topHeight};

        Vec3D axis = top-tip;
        Mdouble length = axis.getLength();
        Vec3D unitAxis = axis/length;

        AxisymmetricIntersectionOfWalls core;
        core.setPosition({0,0,screwTipHeight});
        core.setAxis(unitAxis);
        core.setAngularVelocity({0,0,rotationRate*2*constants::pi});
        core.addObject({-1,0,0},{screwCoreRadius,0,0});
        core.addObject({0,0,-1},{0,0,length});
        core.setSpecies(speciesHandler.getLastObject());
        wallHandler.copyAndAddObject(core);

        //Screw(Vec3D start, Mdouble l, Mdouble r, Mdouble n, Mdouble omega, Mdouble thickness);
        Screw s({0,0,0}, length, screwRadius, screwPitchNum, screwRotationRate, screwThickness);
        s.setPosition(core.getPosition());
        s.setOrientation(core.getOrientation());
        s.setAngularVelocity(core.getAngularVelocity());
        s.setSpecies(core.getSpecies());
        wallHandler.copyAndAddObject(s);
    }

    void addParticlesAtWall() {
        /* Simple run settings
         * Nx*Ny*Nz particles are created evenly spaced between [xmin,xmax]*[ymin,ymax]*[zmin,zmax] and checked for contact with the screw
         */
        const ParticleSpecies* s = speciesHandler.getObject(0);
        particleHandler.clear();

        SphericalParticle p0;
        p0.setSpecies(s);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setRadius(3.0*particleRadius);
        SphericalParticle p1;
        p1.setSpecies(s);
        p1.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p1.setRadius(particleRadius);

        //domain size d
        Vec3D d = Vec3D(getXMax() - getXMin(), getYMax() - getYMin(), getZMax() - getZMin());
        //number of particles that fit in domain
        Mdouble Nx = floor(d.X / (2.0 * particleRadius));
        Mdouble Ny = floor(d.Y / (2.0 * particleRadius));
        Mdouble Nz = floor(d.Z / (2.0 * particleRadius));

        Mdouble distance;
        Vec3D normal;
        Vec3D p;
        Mdouble minDistance;
        int wallIndex;
        unsigned counter = 0;
        for (p.X = getXMin() + particleRadius; p.X < getXMax(); p.X += 2.0 * particleRadius)
            for (p.Y = getYMin() + particleRadius; p.Y < getYMax(); p.Y += 2.0 * particleRadius)
                for (p.Z = getZMin() + particleRadius; p.Z < getZMax(); p.Z += 2.0 * particleRadius)
                {
                    minDistance = p0.getRadius();
                    p0.setPosition(p);
                    for (auto w : wallHandler) {
                        //if touching the wall
                        if (w->getDistanceAndNormal(p0, distance, normal) && distance<minDistance)
                        {
                            minDistance=distance;
                            wallIndex = w->getIndex();
                            if (distance<0) break;
                            p1.setPosition(p0.getPosition()+(distance-particleRadius)*normal);
                        }
                    }
                    if (minDistance<p0.getRadius() && minDistance>0)
                    {
                        p1.setVelocity({0,0,(Mdouble)wallIndex});
                        particleHandler.copyAndAddObject(p1);
                        counter++;
                    }
                }
        std::cout << "Inserted particles: " << counter << std::endl;
    }

    void addParticles() {
        const ParticleSpecies* s = speciesHandler.getObject(0);
        particleHandler.clear();

        SphericalParticle p0;
        p0.setSpecies(s);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setRadius(particleRadius);

        Mdouble distance;
        Vec3D normal;
        Vec3D p;
        Mdouble minDistance;
        for (p.Z = getZMin() + particleRadius; p.Z < 1.1*getZMax(); p.Z += 2.0 * particleRadius*particlePolydispersity) //Changed Cubic to HCP here
            for (p.X = getXMin() + particleRadius; p.X < getXMax(); p.X += 2.0 * particleRadius*particlePolydispersity)
                for (p.Y = getYMin() + particleRadius; p.Y < getYMax(); p.Y += 2.0 * particleRadius*particlePolydispersity)
                {
                    bool touch = false;
                    p0.setPosition(p);
                    for (auto w : wallHandler) {
                        //if touching the wall
                        if (w->getDistanceAndNormal(p0, distance, normal))
                        {
                            touch = true;
                            break;
                        }
                    }
                    if (!touch) {
                        particleHandler.copyAndAddObject(p0);
                        if (particleHandler.getNumberOfObjects() >= maxParticleNumber) {
                            logger(INFO,"Added maximum number of particles (%)",particleHandler.getNumberOfObjects());
                            return;
                        }
                        p0.setRadius(random.getRandomNumber(particleRadius,particleRadius*particlePolydispersity));
                    }
                }
        logger(INFO,"Added % particles",particleHandler.getNumberOfObjects());
    }



};

/**
 * Instantiates the above class and calls solve
 */
int main()
{
    NautaMixer mixer;
    mixer.setName("NautaMixer");
    mixer.setTimeMax(5);
    mixer.solve();
}
