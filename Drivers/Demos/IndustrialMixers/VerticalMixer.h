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

#include "Mercury3D.h"
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Boundaries/CubeInsertionBoundary.h>

/** This code creates a cylindrical container, inserts particles and lets them settle.
*/
struct VerticalMixer : public Mercury3D{

    Mdouble particleRadius_ = 1.5e-3;
    Mdouble drumRadius_ = 30.0*particleRadius_;
    Mdouble drumLength_ = 60.0*particleRadius_;
    Vec3D angularVelocity = Vec3D(1,0,0);
    Mdouble particleNumber_ = 10;
    Mdouble initialParticleOverlap_ = 0;
    //two options for plotting walls
    bool prettyWalls_ = false;
    bool haveOuterWalls = true;

    void setupInitialConditions() override {
        //set general properties
        setMax(Vec3D(0.5*drumLength_, drumRadius_, drumRadius_));
        setMin(-getMax());
        setGravity(Vec3D(0, 0, -9.8));

        //define material/contact properties
        LinearViscoelasticSlidingFrictionSpecies species;
        species.setHandler(&speciesHandler);
        species.setDensity(2000);
        Mdouble mass = species.getMassFromRadius(particleRadius_);
        species.setCollisionTimeAndRestitutionCoefficient(0.005, 0.5, mass);
        species.setSlidingFrictionCoefficient(0.5);
        species.setSlidingStiffness(2.0/7.0*species.getStiffness());
        species.setSlidingDissipation(2.0/7.0*species.getDissipation());
        speciesHandler.copyAndAddObject(species);

        //set timestep
        setTimeStep(species.getCollisionTime(mass)/15);
        //14*collisiontime for nice flow for
        setSaveCount(200);

        //enforce writing wall output
        setWallsWriteVTK(FileType::MULTIPLE_FILES);
        //if (!prettyWalls_)
        //setParticlesWriteVTK(true);

        Mdouble cornerLength = 8.0 * particleRadius_;
        if (!prettyWalls_) {
            //define Walls
            AxisymmetricIntersectionOfWalls outer;
            outer.setSpecies(speciesHandler.getLastObject());
            outer.setAxis(Vec3D(1, 0, 0));
            outer.addObject(Vec3D(1, 0, 0), Vec3D(drumRadius_, 0, 0));
            wallHandler.copyAndAddObject(outer);

            AxisymmetricIntersectionOfWalls leftCorner;
            leftCorner.setSpecies(speciesHandler.getLastObject());
            leftCorner.setAxis(Vec3D(1, 0, 0));
            leftCorner.addObject(Vec3D(1, 0, -1), Vec3D(drumRadius_, 0, -0.5 * drumLength_ + cornerLength));
            wallHandler.copyAndAddObject(leftCorner);

            AxisymmetricIntersectionOfWalls rightCorner;
            rightCorner.setSpecies(speciesHandler.getLastObject());
            rightCorner.setAxis(Vec3D(1, 0, 0));
            rightCorner.addObject(Vec3D(1, 0, 1), Vec3D(drumRadius_, 0, 0.5 * drumLength_ - cornerLength));
            wallHandler.copyAndAddObject(rightCorner);

            InfiniteWall leftEnd;
            leftEnd.setSpecies(speciesHandler.getLastObject());
            leftEnd.set(Vec3D(-1, 0, 0), Vec3D(-0.5 * drumLength_, 0, 0));
            wallHandler.copyAndAddObject(leftEnd);

            InfiniteWall rightEnd;
            rightEnd.setSpecies(speciesHandler.getLastObject());
            rightEnd.set(Vec3D(1, 0, 0), Vec3D(0.5 * drumLength_, 0, 0));
            wallHandler.copyAndAddObject(rightEnd);

            addBlades();

        } else {
            Mdouble thickness = particleRadius_;
            setMax(getMax()+Vec3D(1,1,1)*2*thickness);
            setMin(getMin()-Vec3D(1,1,1)*2*thickness);
            
            if (haveOuterWalls) {
                //define Walls
                AxisymmetricIntersectionOfWalls outer;
                outer.setSpecies(speciesHandler.getLastObject());
                outer.setAxis(Vec3D(1, 0, 0));
                outer.addObject(Vec3D(1, 0, 0), Vec3D(drumRadius_, 0, 0));
                outer.addObject(Vec3D(-1, 0, 0), Vec3D(drumRadius_ + thickness, 0, 0));
                outer.addObject(Vec3D(0, 0, 1), Vec3D(0, 0, -0.5 * drumLength_ + cornerLength));
                outer.addObject(Vec3D(0, 0, -1), Vec3D(0, 0, 0.5 * drumLength_ - cornerLength));
                wallHandler.copyAndAddObject(outer);

                AxisymmetricIntersectionOfWalls leftCorner;
                leftCorner.setSpecies(speciesHandler.getLastObject());
                leftCorner.setAxis(Vec3D(1, 0, 0));
                leftCorner.addObject(Vec3D(1, 0, -1), Vec3D(drumRadius_, 0, -0.5 * drumLength_ + cornerLength));
                leftCorner.addObject(Vec3D(-1, 0, 1),
                                     Vec3D(drumRadius_ + thickness, 0, -0.5 * drumLength_ + cornerLength));
                leftCorner.addObject(Vec3D(0, 0, -1), Vec3D(0, 0, -0.5 * drumLength_ + cornerLength));
                leftCorner.addObject(Vec3D(1, 0, 0), Vec3D(drumRadius_ - cornerLength, 0, 0));
                wallHandler.copyAndAddObject(leftCorner);

                AxisymmetricIntersectionOfWalls rightCorner;
                rightCorner.setSpecies(speciesHandler.getLastObject());
                rightCorner.setAxis(Vec3D(1, 0, 0));
                rightCorner.addObject(Vec3D(1, 0, 1), Vec3D(drumRadius_, 0, 0.5 * drumLength_ - cornerLength));
                rightCorner.addObject(Vec3D(-1, 0, -1),
                                      Vec3D(drumRadius_ + thickness, 0, 0.5 * drumLength_ - cornerLength));
                rightCorner.addObject(Vec3D(0, 0, 1), Vec3D(0, 0, 0.5 * drumLength_ - cornerLength));
                rightCorner.addObject(Vec3D(1, 0, 0), Vec3D(drumRadius_ - cornerLength, 0, 0));
                wallHandler.copyAndAddObject(rightCorner);

                AxisymmetricIntersectionOfWalls leftEnd;
                leftEnd.setSpecies(speciesHandler.getLastObject());
                leftEnd.setAxis(Vec3D(1, 0, 0));
                leftEnd.addObject(Vec3D(0, 0, 1), Vec3D(0, 0, -0.5 * drumLength_ - thickness));
                leftEnd.addObject(Vec3D(0, 0, -1), Vec3D(0, 0, -0.5 * drumLength_));
                leftEnd.addObject(Vec3D(-1, 0, 0), Vec3D(drumRadius_ - cornerLength, 0, 0));
                wallHandler.copyAndAddObject(leftEnd);

                AxisymmetricIntersectionOfWalls rightEnd;
                rightEnd.setSpecies(speciesHandler.getLastObject());
                rightEnd.setAxis(Vec3D(1, 0, 0));
                rightEnd.addObject(Vec3D(0, 0, -1), Vec3D(0, 0, 0.5 * drumLength_ + thickness));
                rightEnd.addObject(Vec3D(0, 0, 1), Vec3D(0, 0, 0.5 * drumLength_));
                rightEnd.addObject(Vec3D(-1, 0, 0), Vec3D(drumRadius_ - cornerLength, 0, 0));
                wallHandler.copyAndAddObject(rightEnd);
            }
            
            addPrettyBlades();
        }

        for (BaseWall* w : wallHandler) {
            w->setAngularVelocity(angularVelocity);
        }

        SphericalParticle p;
        p.setSpecies(speciesHandler.getLastObject());
        CubeInsertionBoundary c; //delete is done in boundaryHandler
        c.set(&p,0,getMin(),getMax(),Vec3D(0,0,0),Vec3D(0,0,0),0.65*particleRadius_-initialParticleOverlap_,1.3*particleRadius_-initialParticleOverlap_);
        c.setHandler(&boundaryHandler);
        while (particleHandler.getSize()<particleNumber_) {
            c.checkBoundaryBeforeTimeStep(this);
            if (particleHandler.getSize()%100==0) logger(INFO,"N %",particleHandler.getSize());
        }
        for (auto p : particleHandler) {
            p->setRadius(p->getRadius()+initialParticleOverlap_);
        }
    }

    void printTime() const override
    {
        logger(INFO,"t %\tN %\tE %\tC %",getTime(),particleHandler.getNumberOfObjects(),
               getKineticEnergy()/getElasticEnergy(),getCentreOfMass().Z);
    }

    virtual void addBlades() {}

    virtual void addPrettyBlades() {}
};

struct VerticalMixerStraightBlades : public VerticalMixer {

    Mdouble bladeWidth_ = 2.0*particleRadius_;
    Mdouble bladeHeight_ = 0.2*drumRadius_;

    void addBlades() {
        IntersectionOfWalls blade;
        blade.setSpecies(speciesHandler.getLastObject());
        blade.createOpenPrism({Vec3D(0,0.5*bladeWidth_,drumRadius_),
                               Vec3D(0,0.5*bladeWidth_,drumRadius_-bladeHeight_),
                               Vec3D(0,-0.5*bladeWidth_,drumRadius_-bladeHeight_),
                               Vec3D(0,-0.5*bladeWidth_,drumRadius_)});
        wallHandler.copyAndAddObject(blade);

        const Vec3D quarterTurn = {2,0,0};
        blade.rotate(quarterTurn);
        wallHandler.copyAndAddObject(blade);
        blade.rotate(quarterTurn);
        wallHandler.copyAndAddObject(blade);
        blade.rotate(quarterTurn);
        wallHandler.copyAndAddObject(blade);
    }
};

struct VerticalMixerAngledBlades : public VerticalMixerStraightBlades {

    Mdouble bladeAngle_ = 0.25*constants::pi;

    void addBlades() {
        Mdouble s = sin(bladeAngle_);
        Mdouble c = cos(bladeAngle_);

        IntersectionOfWalls blade;
        blade.setSpecies(speciesHandler.getLastObject());
        blade.createOpenPrism({Vec3D(0.5*bladeWidth_*s,0.5*bladeWidth_*c,drumRadius_),
                               Vec3D(0.5*bladeWidth_*s,0.5*bladeWidth_*c,drumRadius_-bladeHeight_),
                               Vec3D(-0.5*bladeWidth_*s,-0.5*bladeWidth_*c,drumRadius_-bladeHeight_),
                               Vec3D(-0.5*bladeWidth_*s,-0.5*bladeWidth_*c,drumRadius_)});
        blade.setPosition(Vec3D(0.2*drumLength_,0,0));
        wallHandler.copyAndAddObject(blade);

        const Vec3D quarterTurn = {2,0,0};
        blade.rotate(quarterTurn);
        blade.setPosition(Vec3D(-0.2*drumLength_,0,0));
        wallHandler.copyAndAddObject(blade);
        blade.rotate(quarterTurn);
        blade.setPosition(Vec3D(0.2*drumLength_,0,0));
        wallHandler.copyAndAddObject(blade);
        blade.rotate(quarterTurn);
        blade.setPosition(Vec3D(-0.2*drumLength_,0,0));
        wallHandler.copyAndAddObject(blade);
    }

    void addPrettyBlades() {
        Mdouble s = sin(bladeAngle_);
        Mdouble c = cos(bladeAngle_);

        IntersectionOfWalls blade;
        blade.setSpecies(speciesHandler.getLastObject());
        blade.createPrism({Vec3D(0.5*bladeWidth_*s,0.5*bladeWidth_*c,drumRadius_),
                               Vec3D(0.5*bladeWidth_*s,0.5*bladeWidth_*c,drumRadius_-bladeHeight_),
                               Vec3D(-0.5*bladeWidth_*s,-0.5*bladeWidth_*c,drumRadius_-bladeHeight_),
                               Vec3D(-0.5*bladeWidth_*s,-0.5*bladeWidth_*c,drumRadius_)});
        //restrict to inside
        Mdouble cornerLength = 8.0 * particleRadius_;
        for (Mdouble angle = -0.2*constants::pi; angle<=0.2001*constants::pi; angle+=0.05*constants::pi) {
            Mdouble s = sin(angle);
            Mdouble c = cos(angle);
            //outer.addObject(Vec3D(1, 0, 0), Vec3D(drumRadius_, 0, 0));
            blade.addObject(Vec3D(0, s, -c), Vec3D(0, -s, c)*drumRadius_);
        }
        const Vec3D quarterTurn = {2,0,0};

        IntersectionOfWalls leftBlade = blade;
        for (Mdouble angle = -0.2*constants::pi; angle<=0.2001*constants::pi; angle+=0.05*constants::pi) {
            Mdouble s = sin(angle);
            Mdouble c = cos(angle);
            //leftCorner.addObject(Vec3D(1, 0, -1), Vec3D(drumRadius_, 0, -0.5 * drumLength_ + cornerLength));
            leftBlade.addObject(Vec3D(1, s, -c), Vec3D(0.2*drumLength_-0.5 * drumLength_ + cornerLength,  -s*drumRadius_, c*drumRadius_));
            //blade.addObject(Vec3D(-1, s, -c), Vec3D(-0.2*drumLength_+0.5 * drumLength_ - cornerLength,  -s*drumRadius_, c*drumRadius_));
        }
        leftBlade.setPosition(Vec3D(-0.2*drumLength_,0,0));
        leftBlade.rotate(quarterTurn);
        wallHandler.copyAndAddObject(leftBlade);
        leftBlade.rotate(quarterTurn);
        leftBlade.rotate(quarterTurn);
        wallHandler.copyAndAddObject(leftBlade);

        IntersectionOfWalls rightBlade = blade;
        for (Mdouble angle = -0.2*constants::pi; angle<=0.2001*constants::pi; angle+=0.05*constants::pi) {
            Mdouble s = sin(angle);
            Mdouble c = cos(angle);
            //rightCorner.addObject(...);
            rightBlade.addObject(Vec3D(-1, s, -c), Vec3D(-0.2*drumLength_+0.5 * drumLength_ - cornerLength,  -s*drumRadius_, c*drumRadius_));
        }
        rightBlade.setPosition(Vec3D(0.2*drumLength_,0,0));
        wallHandler.copyAndAddObject(rightBlade);
        rightBlade.rotate(quarterTurn);
        rightBlade.rotate(quarterTurn);
        wallHandler.copyAndAddObject(rightBlade);
    }
};
