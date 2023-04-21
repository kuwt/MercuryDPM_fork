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
#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>

#include <Walls/InfiniteWall.h>
#include "Boundaries/CubeInsertionBoundary.h"
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include "Walls/TriangleWall.h"
#include <Boundaries/CubeDeletionBoundary.h>

using helpers::writeToFile;

#include <map>
#include <random>

using constants::pi;
using mathsFunc::cubic;
enum class RollerType : unsigned char {
    NO_ROLLER = 0,
    FORWARD_ROTATING_ROLLER = 1,
    COUNTER_ROTATING_ROLLER = 2,
    FORWARD_ROTATING_ROLLER_WITH_BLADE = 3
};

class Sintering_S0 : public Mercury3D
{
public:

    Sintering_S0 (ParticleSpecies& particleSpecies, BaseSpecies& wallSpecies)
    {
        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions
        speciesHandler.copyAndAddObject(particleSpecies); //wall-wall interactions
        speciesHandler.getMixedObject(0,1)->mixAll(&wallSpecies,&wallSpecies); //particle-wall interactions (particle-roller)
        //-------------------------------------------------
        //initialTime = getTimeMax();
    }

    void setupInitialConditions() override {

        //++++++++++To store wall information++++++++++++++
        logger(INFO,"Creating a base wall at z=0");
        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0));
        w.set({0,0,-1},{0,0,getZMin()});
        wallHandler.copyAndAddObject(w);
        w.set({0,-1,0},{0,getYMin(),0});
        wallHandler.copyAndAddObject(w);
        w.set({0,1,0},{0,getYMax(),0});
        wallHandler.copyAndAddObject(w);
        w.set({-1,0,0},{getXMin(),0,0});
        wallHandler.copyAndAddObject(w);

        logger(INFO,"Creating a wall at x=Max");
        TriangleWall wall, wall2;
        //auto species = speciesHandler.getObject(0);
        wall.setSpecies(speciesHandler.getObject(0));
        wall.setVertices(Vec3D(getXMax(),getYMax(),0.0),Vec3D(getXMax(),getYMax(),getZMax()),Vec3D(getXMax(),getYMin(),getZMax()));

        wall2.setSpecies(speciesHandler.getObject(0));
        wall2.setVertices(Vec3D(getXMax(),getYMin(),getZMax()),Vec3D(getXMax(),getYMin(),0),Vec3D(getXMax(),getYMax(),0));

        wallHandler.copyAndAddObject(wall);
        wallHandler.copyAndAddObject(wall2);

        //++++++++++To store particle information++++++++++++++
        logger(INFO,"Adding particles ...");
        /* Introduce the InsertionBoundary */
        ThermalParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));

        insb = new CubeInsertionBoundary;

        insb->set(
                insertionBoundaryParticle,
                10,
                Vec3D(getXMin(), getYMin(), getZMin()),
                Vec3D(0.9*getXMax(), 0.9*getYMax(), getZMax()*3.0),
                Vec3D(0, 0, 0),
                Vec3D(0, 0, 0)
        );
        PSD psd;
        psd.setDistributionUniform(radius*0.98, radius*1.1, 50);
        insb->setPSD(psd);
        insb = boundaryHandler.copyAndAddObject(insb);

        //++++++++++To delete particles at specific position++++++++++++++
        delb = new DeletionBoundary;
        delb->set(Vec3D(0,0,1), getZMax());
        delb = boundaryHandler.copyAndAddObject(delb);
        delb->set(Vec3D(1,0,0), getXMax());
        delb = boundaryHandler.copyAndAddObject(delb);

        logger(INFO,"Inserted % particles",boundaryHandler.getNumberOfObjects());
    }

    //remove particles according to (max.) outflow rate
    void actionsAfterTimeStep() override
    {
        //store time when rolling should start
        static Mdouble rollerTime = 0;
        //make decisions based on compression stage
        static unsigned compressionStage_ = 0;

        if (compressionStage_==0)
        {
            // Functions after this statement only get executed every 100th time step (if counter==100)
            static unsigned counter = 0;
            if (++counter != 100) return;
            else counter = 0;

            //Mdouble timeToRoll = getTime();

            //if (getKineticEnergy() < 5e-3 * getElasticEnergy())
            if (getTime() > 0.6)
            {
                //change to compression state
                compressionStage_ = 1;

                //set roller time
                rollerTime = getTime() + 2.0*getZMax() / rollerVelocity_;

                //set final time such that the roller completes the domain
                setTimeMax((rollerTime + (getXMax() - getXMin() - rollerRadius_) / rollerVelocity_) + 0.25);

                logger(INFO, "Starting compression at t=%, rollerTime=%, tMax=%", getTime(), rollerTime, getTimeMax());

                //create roller
                roller = new AxisymmetricIntersectionOfWalls;
                roller->setSpecies(speciesHandler.getObject(1));
                roller->setPosition({rollerRadius_, 0, rollerRadius_ + 2.0*getZMax() + compaction_});
                roller->setAxis({0, 1, 0});//setAxis
                roller->addObject({-1, 0, 0}, {rollerRadius_, 0, 0});
                roller->setVelocity({0,0,-rollerVelocity_});//first move the roller down
                roller = wallHandler.copyAndAddObject(roller);
            }
        } else {
            if (getTime() > rollerTime ) {
                logger(INFO, "Start forward movement");
                // set velocity
                wallHandler.getObject(1)->setVelocity({rollerVelocity_,0,0});
                wallHandler.getLastObject()->setVelocity({rollerVelocity_,0,0});
                //set angular velocity
                if (rollerType_ == RollerType::COUNTER_ROTATING_ROLLER)
                {
                    wallHandler.getObject(1)->setAngularVelocity({0,-rollerVelocity_/rollerRadius_,0});
                } else {
                    wallHandler.getObject(1)->setAngularVelocity({0,rollerVelocity_ /rollerRadius_,0});
                }
                logger(INFO, "a=%",wallHandler.getObject(1)->getAngularVelocity());
                logger(INFO, "z=%",wallHandler.getObject(1)->getPosition().Z-rollerRadius_);
                logger(INFO, "z=% %",compaction_,getZMax());

                // make sure this is set only once
                rollerTime = getTimeMax();
                boundaryHandler.removeObject(insb->getIndex());

                TrackingRoller = rollerTime;
            }
        }
    }

    void addRoller (RollerType rollerType, Mdouble rollerRadius, Mdouble rollerVelocity, Mdouble compaction, Mdouble bladeThickness=50e-6, Mdouble bladeRollerDistance=50e-6)
    {
        rollerType_ = rollerType;
        rollerRadius_ = rollerRadius;
        rollerVelocity_ = rollerVelocity;
        compaction_ = compaction;
        bladeThickness_ = bladeThickness;
        bladeRollerDistance_ = bladeRollerDistance;
    }

    void printTime() const override
    {
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << " tmax " << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << "  EneRatio " << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
                  << "  KineticEnergy " << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()
                  << "  ElasticEnergy " << std::setprecision(3) << std::left << std::setw(6) << getElasticEnergy()
                  << "  TrackingRAdholler " << std::setprecision(3) << std::left << std::setw(6) << TrackingRoller
                  << std::endl;
    }

    Mdouble radius = 0.5;

    private:

    Mdouble relativeBasalLayerThickness_;
    Mdouble rollerRadius_;
    Mdouble rollerVelocity_;
    Mdouble compaction_;
    Mdouble bladeThickness_;
    Mdouble bladeRollerDistance_;
    RollerType rollerType_ = RollerType::NO_ROLLER;
    CubeInsertionBoundary* insb;
    AxisymmetricIntersectionOfWalls* roller;

    DeletionBoundary* delb;
    Mdouble TrackingRoller;

    unsigned int TimeRoller;

    unsigned int num_inserted;
    double mass_inserted, vol_inserted;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{

    //define domain size
    Mdouble domainLength = 0.15; //[mm]
    Mdouble domainWidth = 0.15; //[mm]
    Mdouble domainDepth = 0.015; //[mm]

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 1.04; //[mg/mm^3]

    Mdouble meanRadius = 0.002; //[mm]//[

    Mdouble mass = density*(4.0/3.0*pi)*cubic(meanRadius);

    //**********************************************************************
    //Mdouble collisionTime = 2.0e-4; //-> I want a stiffness= 5e6 [kg/s^2]
    Mdouble stiffness = 0.03; //[mg/mu s^2]
    Mdouble gamma = 1.0e-6; //[mg/mu s]
    //Mdouble restitution = 0.2;
    //**********************************************************************

    //----------------------------------------------
    Mdouble friction = 0.55;
    Mdouble phi_f =0.05; //Penetration Depth [I need it!] 5%

    //----------------------------------------------
    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //Return Stiffness and dissipation
    //particleSpecies.setCollisionTimeAndRestitutionCoefficient(collisionTime, restitution, mass);
    //particleSpecies.setStiffnessAndRestitutionCoefficient(stiffness,restitution,mass);

    Mdouble k2 = stiffness;//LoadingStiffness
    Mdouble k1 = 5.0*k2; //Unloading stiffness
    Mdouble kc = 0.5*k2; //Cohesive stiffness

    //-------->[Set Plastic Contribution]
    particleSpecies.setPlasticParameters(k2, k1, kc, phi_f);

    //----------------------------------------------
    particleSpecies.setDissipation(gamma);
    particleSpecies.setDensity(density);
    //----------------------------------------------
    Mdouble k_s = 0.2*particleSpecies.getLoadingStiffness();
    Mdouble k_r = 0.1*particleSpecies.getLoadingStiffness();
    Mdouble k_o = 0.1*particleSpecies.getLoadingStiffness();

    Mdouble miu_s = 1.0*friction;
    Mdouble miu_r = 0.1*friction;
    Mdouble miu_o = 0.1*friction;

    Mdouble gamma_s = 0.2*gamma;
    Mdouble gamma_r= 0.05*gamma;
    Mdouble gamma_o= 0.05*gamma;

    //----------------------------------------------
    //-------->[Set Tangential Contribution]
    particleSpecies.setSlidingStiffness(k_s);
    particleSpecies.setRollingStiffness(k_r);
    particleSpecies.setTorsionStiffness(k_o);

    particleSpecies.setSlidingFrictionCoefficient(miu_s);
    particleSpecies.setRollingFrictionCoefficient(miu_r);
    particleSpecies.setTorsionFrictionCoefficient(miu_o);

    particleSpecies.setSlidingDissipation(gamma_s);
    particleSpecies.setRollingDissipation(gamma_r);
    particleSpecies.setTorsionDissipation(gamma_o);

    //----------------------------------------------
    //-------->[Set Adhesive Contribution]
    //Set Adhesive properties:
    Mdouble K_adh = k2;//0.5*k2;//1.0;
    Mdouble f_adh_max = 1e-10*k2;//5.0*(meanRadius*2.0*meanRadius*2.0)*5.0; ///set adhesive force max based on Bond number [Hao]

    particleSpecies.setAdhesionStiffness(K_adh);
    particleSpecies.setAdhesionForceMax(f_adh_max);

    //----------------------------------------------
    //-------->[Set Sinter Contribution]
    //Set sinter parameters:
    //particleSpecies.setSinterType(SINTERTYPE::CONSTANT_RATE);
    //particleSpecies.setSinterRate(0.0);
    //particleSpecies.setSinterAdhesion(1200e-9/radius);

    //----------------------------------------------
    //create the species for particle-wall interactions
    ThermalSinterLinFrictionReversibleAdhesiveMixedSpecies wallSpecies = particleSpecies;

    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    Sintering_S0 pb (particleSpecies, wallSpecies);

    pb.radius = meanRadius;
    logger(INFO,"Gravitational compression % of radius per layer",100*mass*9.8/particleSpecies.getLoadingStiffness()/25e-6);

    //----------------------------------------------
    pb.setGravity({0,0,-9.81}); //mm/mu s^2 [Take care]

    pb.setMin({0,-0.5*domainWidth,0});
    pb.setMax({domainLength,0.5*domainWidth,domainDepth});
    //----------------------------------------------

    //----------------------------------------------
    Mdouble tc1 = particleSpecies.getCollisionTime(mass)/50.0;

    Mdouble tc2 = particleSpecies.computeTimeStep(mass)*2.0;

    //----------------------------------------------
    Mdouble GravityTime = sqrt((meanRadius*2.0)/9.81e-9);
    Mdouble CollistionTime = particleSpecies.getCollisionTime(mass);
    Mdouble TimeStepe = particleSpecies.getCollisionTime(mass)/50;

    //----------------------------------------------
    pb.setTimeStep(TimeStepe);
    pb.setTimeMax(2.0);

    //----------------------ROLLER---------------------------
    Mdouble rollerRadius = 0.035; //[mm]
    Mdouble rollerVelocity = 0.08; // mm/mu s
    Mdouble compaction = 1.8*domainDepth;
    pb.addRoller(RollerType::COUNTER_ROTATING_ROLLER,rollerRadius,rollerVelocity,compaction);

    //----------------------------------------------
    pb.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
    pb.setParticlesWriteVTK(true);

    //pb.setSaveCount(100);

    //pb.setNumberOfDomains({1,1,2});
    pb.setName("S0_InitialConditions");
    pb.setXBallsAdditionalArguments("-solidf -v0");
    pb.setFileType(FileType::ONE_FILE);


    pb.solve();

    return 0;
}