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

using constants::pi;
using mathsFunc::cubic;

class LaserBeam : public Mercury3D
{
private:

    Mdouble radius_ = 0.0;
    Mdouble deltaR_ = 0.0; //[check] Thermal expansion coefficient
    Mdouble numLayers_ = 1.0; //1.0; By default at least one layer.

    Mdouble startingTemp_ = 0.0; //[C]
    Mdouble gradientTemp_ = 0.0; //[C/s] Heating rate
    Mdouble meltTemp_ = 0.0; //[C]
    Mdouble maxTemp_ = 0.0;

    Mdouble deltaTemperature = 0.0;
    Mdouble timeCooling_ = 0.0; //[s]

    Vec3D setLaserPosition_;

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint_;
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* wallSpecies; //pointer to the wall properties

public:

    LaserBeam (Mdouble radius, Mdouble deltaR, Mdouble numLayers, ParticleSpecies& particleSpecies, Mdouble  startingTemp, Mdouble gradientTemp,
               Mdouble meltTemp, Mdouble maxTemp, Mdouble timeCooling)
    {
        radius_ = radius;
        deltaR_ = deltaR;
        numLayers_ = numLayers;

        startingTemp_ =  startingTemp;
        gradientTemp_ = gradientTemp;
        meltTemp_ = meltTemp;
        maxTemp_ = maxTemp;
        timeCooling_ = timeCooling;

        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions
        //-------------------------------------------------

        //WallSpecies. This parameters are always constant.
        wallSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        wallSpecies-> setDensity(4000);
//        wallSpecies->setDissipation(particleSpecies->getDissipation());
        const Mdouble YoungM = 1.65e9; //[Pa] Young's Modulus
        Mdouble k_wall = YoungM;
        wallSpecies->setPlasticParameters(100*k_wall,200*k_wall, 0.0, 0.001);

        //wallParticleSpecies = speciesHandler.getMixedObject(particleSpecies,wallSpecies)->mixAll(particleSpecies,wallSpecies);
        speciesHandler.getMixedObject(0,1)->mixAll(&particleSpecies,wallSpecies); //particle-wall interactions (particle-roller)

        particleSpecies_ = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
        particleSpeciesAtMeltingPoint_ = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies_);
    }

    void setupInitialConditions() override {

        particleSpeciesAtMeltingPoint_->setLoadingStiffness(particleSpecies_->getLoadingStiffness());
        particleSpeciesAtMeltingPoint_->setSinterAdhesion(particleSpecies_->getSinterAdhesion());

        //++++++++++To store particle information++++++++++++++
        logger(INFO,"Adding particles ...");
        /* Introduce the InsertionBoundary */
        ThermalParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));

        //--------------------------------------------------
        hGridRebuild();
        ThermalParticle P;
        P.setRadius(radius_);
        P.setSpecies(speciesHandler.getObject(0));

//        std::array<Vec3D, 4> layerPos = {Vec3D(1, 1, 0), Vec3D(1, -1, 0), Vec3D(-1, -1, 0), Vec3D(-1, 1, 0)};
        std::array<Vec3D, 36> layerPos = {Vec3D(1, 1, 0), Vec3D(1, -1, 0), Vec3D(1, -3, -0), Vec3D(-1, -3, 0),
                                          Vec3D(1, -5, 0), Vec3D(-1, -5, 0),
                                          Vec3D(-1, -1, 0), Vec3D(-1, 1, 0),
                                          Vec3D(3, 1, 0), Vec3D(3, -1, 0),
                                          Vec3D(3, -3, 0), Vec3D(3, -5, 0),
                                          Vec3D(5, 1, 0), Vec3D(5, -1, -0),
                                          Vec3D(5, -3, 0), Vec3D(5, -5, 0),
                                          Vec3D(7, 1, 0), Vec3D(7, -1, 0),
                                          Vec3D(7, -3, 0), Vec3D(7, -5, 0),
                                          Vec3D(7, -7, 0), Vec3D(5, -7, 0),
                                          Vec3D(3, -7, 0), Vec3D(1, -7, 0),
                                          Vec3D(-1, -7, 0),Vec3D(9, 1, 0),
                                          Vec3D(9, -1, 0),Vec3D(9, -3, 0),
                                          Vec3D(9, -5, 0),Vec3D(9, -7, 0),
                                          Vec3D(9, -9, 0),
                                          Vec3D(7, -9, 0),Vec3D(5, -9, 0),
                                          Vec3D(3, -9, 0),Vec3D(1, -9, 0),
                                          Vec3D(-1, -9, 0)};
//        double posFirst = 11*radius_;
        double posFirst = radius_;
//        double posFirst = 1.0;
        for (unsigned i = 0; i < numLayers_; i++) {
            Vec3D center = Vec3D(getXMin()+2*radius_,
                                 getYMin()+2*radius_,
                                 getZMin()+ i * 2 * radius_ + posFirst);
            for (auto pos : layerPos) {
                P.setPosition(center + radius_ * pos*1.0);
                P.setTemperature(startingTemp_);
                P.setVelocity(Vec3D(0.0, 0.0, 0.0));
                particleHandler.copyAndAddObject(P);
                /// fix the last particles in space
//                if (i == nLayers - 1) {
//                    particleHandler.getLastObject()->fixParticle();
//                    /// give an initial kick to the last layer of particles
////                    particleHandler.getLastObject()->setPosition(P.getPosition() + (Vec3D(0.0, -0.05 * radius_, 0.0)));
//                }
            }
        }

        std::cout << std::endl;
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
        std::cout << std::endl;

        InfiniteWall w;
        w.setSpecies(wallSpecies);
        w.set({0,0,-1},{0,0,getZMin()});
        wallHandler.copyAndAddObject(w);
//        w.set({0,0,1},{0,0,getZMax()});
//        wallHandler.copyAndAddObject(w);
        w.set({0,-1,0},{0,getYMin(),0});
        wallHandler.copyAndAddObject(w);
//        w.set({0,1,0},{0,getYMax(),0});
//        wallHandler.copyAndAddObject(w);
        w.set({-1,0,0},{getXMin(),0,0});
        wallHandler.copyAndAddObject(w);
//        w.set({1,0,0},{getXMax(),0,0});
//        wallHandler.copyAndAddObject(w);
    }

    //--------------------------------------------------
    Mdouble getTemperatureProfile()const
    {
        Mdouble heatingTemperature = startingTemp_ + gradientTemp_*getTime();
        Mdouble coolingTemperature = maxTemp_;

        if(getTime()>=timeCooling_){
            coolingTemperature = startingTemp_ + gradientTemp_*(getTimeMax() - getTime())*0.7;
        }
        return std::min(std::min(heatingTemperature,coolingTemperature),maxTemp_);
    }

    Mdouble getSurfaceHeight() const
    {
        double height = 0;

        for (const BaseParticle* p : particleHandler)
        {
            double newHeight = p->getPosition().Z;
            if (height<newHeight) height = newHeight;
        }
        return height+particleHandler.getMeanRadius();
    }

    //Change particle temperature according to temperature profile
    void setTemperature(Mdouble tempStep)
    {
        static Mdouble oldTemperature = startingTemp_;
        deltaTemperature = tempStep-oldTemperature;
        Mdouble deltaCooling = tempStep - meltTemp_;

        Mdouble factorRadius = 1.0 - (deltaR_ * deltaTemperature);

        //change species properties
        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltTemp_-tempStep)/gradientTemp_));
        Mdouble oldLoadingStiffness = particleSpecies_->getLoadingStiffness();


//        for (BaseParticle* p : particleHandler)
//        {
//            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
//            tp->setTemperature(tempStep);
//        }
//
//        for (BaseParticle* p : particleHandler)
//        {
//            if (p->getSpecies()==particleSpecies_)
//            {
//                p->setRadius((factorRadius * p->getRadius()));
//            }
//        }
        //add temperature above a certain height
        Mdouble lb_z = getSurfaceHeight() - 2.0 * particleHandler.getMeanRadius();
        Mdouble lb_x = setLaserPosition_.X;
        Mdouble lb_y = setLaserPosition_.Y;

        //heat particles in top layer - Heat Transfer
        for (BaseParticle *p : particleHandler) {
            if (p->getPosition().Z > lb_z)
            {
                auto *tp = dynamic_cast<ThermalParticle *>(p);

//                particleSpecies_->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint_->getLoadingStiffness());
//                particleSpecies_->setSinterAdhesion(1000.0*particleSpeciesAtMeltingPoint_->getLoadingStiffness());
//                particleSpecies_->setFluidity(6.0);
//                particleSpecies_->setSeparationDis(7.0e-7);

//                logger(INFO, "Index= %", p->getIndex());
//                //logger(INFO, "Stiffness= %", particleSpecies_->getLoadingStiffness());
                //
                tp->setTemperature(tempStep); //Change the temperature of the particle
//                p->setRadius((factorRadius * p->getRadius())); //change the radius of the heated particle
//                particleHandler.getObject(p->getIndex())->setSpecies(particleSpecies_);
            } else{

//                particleSpecies_->setLoadingStiffness(particleSpeciesAtMeltingPoint_->getLoadingStiffness());
//                particleSpecies_->setSinterAdhesion(0.000001*particleSpeciesAtMeltingPoint_->getLoadingStiffness());
//                particleSpecies_->setFluidity(0.0);
//                particleSpecies_->setSeparationDis(0.0);
//                particleHandler.getObject(p->getIndex())->setSpecies(particleSpecies_);
            }
        }

//        particleSpecies_->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint_->getLoadingStiffness());

        //for decreasing temperature, change the maxOverlap
//        if (deltaCooling<0.0)
//        {
//            for (BaseInteraction* cBase : interactionHandler)
//            {
//                auto c =  dynamic_cast<SinterLinInteraction*>(cBase);
//                Mdouble unloadingStiffness = c->getUnloadingStiffness();
//                c->setMaxOverlap(c->getMaxOverlap()
//                                 *(unloadingStiffness-oldLoadingStiffness)
//                                 /(unloadingStiffness-particleSpecies_->getLoadingStiffness())
//                );
//            }
//        }
        oldTemperature = tempStep;
    }

    void setLaserPosition()
    {
        //This vector will change eventually.
        Vec3D getCenter;
        getCenter.X = (getXMax()-getXMin())/2.0;
        getCenter.Y = (getYMax()-getYMin())/2.0;
        getCenter.Z = (getZMax()-getZMin())/2.0;

        setLaserPosition_.X =  getCenter.X;
        setLaserPosition_.Y = getCenter.Y;
        setLaserPosition_.Z = getSurfaceHeight();
    }

    //--------------------------------------------------
    void actionsAfterTimeStep() override
    {
        setLaserPosition();
        setTemperature(getTemperatureProfile()); //With Temperature profile. i.e., dilatometer

    }

    //--------------------------------------------------
    double getInfo(const BaseParticle& p) const override
    {
        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-startingTemp_)/(maxTemp_-startingTemp_);
    }

    //System
    Mdouble getThermalEnergy() const //ToDo: Check units.
    {
        Mdouble eneEnergy;

        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;

            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            eneEnergy += (tp->getTemperature()-startingTemp_) * tp->getMass() * particleSpecies_->getHeatCapacity();
//            eneEnergy += deltaTemperature * tp->getMass() * particleSpecies->getHeatCapacity();
        }
        return eneEnergy;
    }

    //--------------------------------------------------
    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
        return sqrt(meanOverlap/meanRadius);
    }

    //    Function to override the output file with specific parameters
    void writeFstatHeader(std::ostream &os) const override
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble Force = abs(dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getForce().X);
        Mdouble meanRadius = particleHandler.getMeanRadius();
        Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;
        Mdouble stress = Force/meanContactRadius;

        Mdouble X = getMeanRelativeContactRadius()/2;
        Mdouble rho_p = 1.0/(1 + std::pow(1 - X,3.0));

        os << getTime() //1
           << " " << getTemperatureProfile() //2
           << " " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getRadius()//3
           << " " << particleSpecies_->getDensity()//4
           << " " << Force//5
           << " " << stress //6
           << " " << particleHandler.getMass()//7
           << " " << particleSpecies_->getLoadingStiffness()//8
           << " " << meanOverlap//9
           << " " << meanContactRadius //10
           << " " <<  particleHandler.getMeanRadius() //11
           << " " <<  getMeanRelativeContactRadius()//12 // Neckgrowth
           << " " << getThermalEnergy()/2//13 [J]
           << " " << getThermalEnergy()/(getTime()) //14 [W/s]
           << " " << getKineticEnergy() //15
           << " " << getElasticEnergy() //16
           << " " << rho_p //17
           << std::endl;
    }

    void printTime() const override
    {
        Mdouble volSystem = (getXMax()-getXMin())* (getYMax()-getYMin())*(getZMax()-getZMin());
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax " <<  std::setprecision(3) << std::left << std::setw(3) << getTimeMax()
                  << ", Temperature " << std::setprecision(3) << std::left << std::setw(6) << getTemperatureProfile()
                  << ", pRadius " << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getMeanRadius()
                  << ", K1: " << std::setprecision(3) << std::left << std::setw(3) << particleSpecies_->getLoadingStiffness()
                  << ", Heat Energy " << std::setprecision(3) << std::left << std::setw(6) << getThermalEnergy()
                  << ", square(delta/R) " << std::setprecision(3) << std::left << std::setw(6) << getMeanRelativeContactRadius()
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "SetSpecificParticlesToHeat";

    Mdouble meanRadius = 3.0e-5; //[mm]//
    Mdouble deltaR = 1.0e-3; //Relate to thermal expasion coefficient
    Mdouble numLayers = 3;
    //define domain size
    Mdouble domainLength = 6*(2*meanRadius); //[mm]
    Mdouble domainWidth = 6*(2*meanRadius); //[mm]
    Mdouble domainDepth = numLayers*(2*meanRadius); //[mm]

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 930.0; //[mg/mm^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 1.650e9; //Young's modulus
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble pDepth = 1.45; //Penetration depth
    const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);

    const Mdouble k1 =0.01*E_reduced*r_effective;
    const Mdouble k2 = 10.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    //----------------------------------------------
    Mdouble setDeltaC = 5.1e-07 ;
    Mdouble setC1 = 25.2;
    Mdouble sinterAdhesion = 0.00001*k1;

    //Thermal variables:
    //PS
    Mdouble setHCapacity = 1320.0; //[J/(kgK)]
    Mdouble setTConductivity = 0.167;  //[W/(mK)]

    Mdouble startingTemp_ = 20.0;
    Mdouble gradientTemp = 440.0;
    Mdouble meltTemp = 40.0;
    Mdouble maxTemp = 70.0;
    Mdouble timeCooling = 0.18;

    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //-------->[Set Plastic Contribution]
    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
    particleSpecies.setDensity(density);
    particleSpecies.setSinterAdhesion(sinterAdhesion);

    //-------------------
    //Sinter parameters
    particleSpecies.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
    particleSpecies.setComplianceZero(1/(2*E_Modulus)); //Book: Adhesive Particles
    particleSpecies.setSurfTension(0.035); // [N/m] or [J/m^2] surface energy density or surface tension.

    //Parameters to calibrate:
    particleSpecies.setFluidity(setC1);
    particleSpecies.setSeparationDis(setDeltaC);

    particleSpecies.setThermalConductivity(setTConductivity);
    particleSpecies.setHeatCapacity(setHCapacity);

    //----------------------------------------------
    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    LaserBeam pb (meanRadius,deltaR,numLayers,particleSpecies,startingTemp_,gradientTemp,
                  meltTemp, maxTemp, timeCooling);

    //----------------------------------------------
    pb.setGravity({0,0,0.0}); //mm/mu s^2 [Take care]
    pb.setMin({0.0,0.0,0.0});
    pb.setMax({domainLength,domainWidth,domainDepth});
    //----------------------------------------------

    Mdouble TimeStep = particleSpecies.getCollisionTime(mass)/50;

    //----------------------------------------------
    pb.setSaveCount(1000);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(1.0);

    logger(INFO,"Time step: %", pb.getTimeStep());
    //----------------------------------------------
//    pb.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
    pb.setParticlesWriteVTK(true);

    pb.setName(setFilename);
    pb.setXBallsAdditionalArguments("-solidf -v0");
    pb.setFileType(FileType::ONE_FILE);

    pb.removeOldFiles();
    pb.setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");
//    pb.setXBallsAdditionalArguments("-solidf -v0 ");

    pb.solve();

    helpers::writeToFile(setFilename + ".gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "plot 'SetSpecificParticlesToHeat.fstat' u 1:12 title 'DEM simulation'  with lines linestyle 2\n"
    );

    helpers::writeToFile(setFilename + "2.gnu",
                         "set xlabel 'a(t)/R'\n"
                         "set ylabel 'Q [J]'\n"
                         "plot 'SetSpecificParticlesToHeat.fstat' u 12:13 title 'Heat energy'  with lines linestyle 2"
    );

    return 0;
}