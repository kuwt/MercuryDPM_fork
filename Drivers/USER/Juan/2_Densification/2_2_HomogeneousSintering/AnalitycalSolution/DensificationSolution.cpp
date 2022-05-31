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

class Densification : public Mercury3D
{
private:

    Mdouble radius_ = 0.0;
    Mdouble deltaR_ = 0.0; //[check] Thermal expansion coefficient
    Mdouble numLayers_ = 0.0;

    Mdouble startingTemp_ = 0.0; //[C]
    Mdouble gradientTemp_ = 0.0; //[C/s] Heating rate
    Mdouble meltTemp_ = 0.0; //[C]
    Mdouble maxTemp_ = 0.0;

    Mdouble deltaTemperature = 0.0;

    Mdouble timeCooling_ = 0.0; //[s]

    Mdouble scalarNormalForce = 0.0;
    Mdouble meanCoordinationNumber = 0.0;

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint_;

public:

    Densification (Mdouble radius, Mdouble deltaR, Mdouble numLayers, ParticleSpecies& particleSpecies, Mdouble  startingTemp, Mdouble gradientTemp,
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
        particleSpecies_ = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
        particleSpeciesAtMeltingPoint_ = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies_);
    }

    void setupInitialConditions() override {

        particleSpeciesAtMeltingPoint_->setLoadingStiffness(particleSpecies_->getLoadingStiffness());

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

        unsigned nLayers = numLayers_;
//        std::array<Vec3D, 9> layerPos = {Vec3D(1, 1, 0), Vec3D(1, -1, 0), Vec3D(-1, -1, 0), Vec3D(-1, 1, 0),
//                                         Vec3D(3, 1, 0), Vec3D(3, -1, 0), Vec3D(3, -3, 0), Vec3D(1, -3, -0), Vec3D(-1, -3, 0)};
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

//        std::array<Vec3D, 36> layerPos = {Vec3D(1, 0, 1), Vec3D(1, 0, -1), Vec3D(1, 0, -3), Vec3D(-1, 0, -3),
//                                         Vec3D(1, 0, -5), Vec3D(-1, 0, -5),
//                                         Vec3D(-1, 0, -1), Vec3D(-1, 0, 1),
//                                         Vec3D(3, 0, 1), Vec3D(3, 0, -1),
//                                          Vec3D(3, 0, -3), Vec3D(3, 0, -5),
//                                          Vec3D(5, 0, 1), Vec3D(5, 0, -1),
//                                          Vec3D(5, 0, -3), Vec3D(5, 0, -5),
//                                          Vec3D(7, 0, 1), Vec3D(7, 0, -1),
//                                          Vec3D(7, 0, -3), Vec3D(7, 0, -5),
//                                          Vec3D(7, 0, -7), Vec3D(5, 0, -7),
//                                          Vec3D(3, 0, -7), Vec3D(1, 0, -7),
//                                          Vec3D(-1, 0, -7),Vec3D(9, 0, 1),
//                                          Vec3D(9, 0, -1),Vec3D(9, 0, -3),
//                                          Vec3D(9, 0, -5),Vec3D(9, 0, -7),
//                                          Vec3D(9, 0, -9),
//                                          Vec3D(7, 0, -9),Vec3D(5, 0, -9),
//                                          Vec3D(3, 0, -9),Vec3D(1, 0, -9),
//                                          Vec3D(-1, 0, -9)};

        double posFirst = radius_;
//        double posFirst = 1.0;
        for (unsigned i = 0; i < numLayers_; i++) {
            Vec3D center = Vec3D(getXMin()+ 2.0*radius_,
                                 getYMin()+ 2.0*radius_,
                                 getZMin()+ i * 1.999 * radius_ + posFirst);
            for (auto pos : layerPos) {
                P.setPosition(center + radius_ * pos*0.999);
//                P.fixParticle();
                P.setTemperature(startingTemp_);
                P.setVelocity(Vec3D(0.0, 0.0, 0.0));
                particleHandler.copyAndAddObject(P);
            }
        }

        std::cout << std::endl;
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
        std::cout << std::endl;
    }

    //--------------------------------------------------
    Mdouble getTemperature()const
    {
        Mdouble heatingTemperature = startingTemp_ + gradientTemp_*getTime();
        Mdouble coolingTemperature = maxTemp_;

        if(getTime()>=getTimeMax()-timeCooling_){
            //Cooling rate is 70%  less than heating rate.
            coolingTemperature = startingTemp_ + (gradientTemp_*0.7)*(getTimeMax() - getTime())*0.7;
        }

        return std::min(std::min(heatingTemperature,coolingTemperature),maxTemp_);
    }

    Mdouble getNewVolume() const
    {
        double heightMax = 0.0;
        Mdouble heightMin = getZMax();

        for (const BaseParticle* p : particleHandler)
        {
            double newHeightMax = p->getPosition().Z;
            double newHeightMin = p->getPosition().Z;

            if (heightMax<newHeightMax) heightMax = newHeightMax;
            if (newHeightMin< heightMin) heightMin = newHeightMin;

        }

        heightMax += particleHandler.getMeanRadius();
        heightMin -= particleHandler.getMeanRadius();

        //+++
        double widthMax = 0.0;
        double widthMin = getXMax();

        for (const BaseParticle* p : particleHandler)
        {
            double newWidthMax = p->getPosition().X;
            double newWidthMin = p->getPosition().X;
            if (widthMax<newWidthMax) widthMax = newWidthMax;
            if (newWidthMin<widthMin) widthMin = newWidthMin;
        }

        widthMax += particleHandler.getMeanRadius();
        widthMin -= particleHandler.getMeanRadius();
        //+++

        double lengthMax = 0.0;
        double lengthMin = getYMax();

        for (const BaseParticle* p : particleHandler)
        {
            double newLengthMax = p->getPosition().Y;
            double newLengthMin = p->getPosition().Y;

            if (lengthMax<newLengthMax) lengthMax = newLengthMax;
            if (newLengthMin<lengthMin)lengthMin = newLengthMin;

        }

        lengthMax += particleHandler.getMeanRadius();
        lengthMin -= particleHandler.getMeanRadius();
        //+++

        return (heightMax-heightMin)*(widthMax-widthMin)*(lengthMax-lengthMin);

    }

    //Change particle temperature according to temperature profile
    void setTemperature(Mdouble tempStep)
    {
        for (BaseParticle* p : particleHandler)
        {
            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            tp->setTemperature(tempStep);
        }

        static Mdouble oldTemperature = startingTemp_;
        deltaTemperature = tempStep-oldTemperature;
        Mdouble deltaCooling = tempStep - meltTemp_;

        Mdouble factorRadius = 1.0 + (deltaR_ * deltaTemperature);

        for (BaseParticle* p : particleHandler)
        {
            if (p->getSpecies()==particleSpecies_)
            {
                p->setRadius((factorRadius * p->getRadius()));
            }
        }
        //change species properties
        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltTemp_-tempStep)/gradientTemp_));
        Mdouble oldLoadingStiffness = particleSpecies_->getLoadingStiffness();
        particleSpecies_->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint_->getLoadingStiffness());

        if(tempStep ==  meltTemp_)
        {
//            particleSpecies_->setLoadingStiffness(0.08*stiffnessFactor*particleSpeciesAtMeltingPoint->getLoadingStiffness()); //Todo!
            particleSpecies_->setSinterAdhesion(0.0009*particleSpeciesAtMeltingPoint_->getLoadingStiffness());

        }else if(tempStep >=  0.9*meltTemp_){

//            particleSpecies_->setLoadingStiffness(0.1*stiffnessFactor*particleSpeciesAtMeltingPoint->getLoadingStiffness()); //Todo!
            particleSpecies_->setSinterAdhesion(0.0005*particleSpeciesAtMeltingPoint_->getLoadingStiffness());

        }
        else{
//            particleSpecies_->setLoadingStiffness(0.3*stiffnessFactor*particleSpeciesAtMeltingPoint->getLoadingStiffness()); //Todo!
            particleSpecies_->setSinterAdhesion(0.0001*particleSpeciesAtMeltingPoint_->getLoadingStiffness());

        }

        //for decreasing temperature, change the maxOverlap
        if (deltaCooling<0.0)
        {
            for (BaseInteraction* cBase : interactionHandler)
            {
                auto c =  dynamic_cast<SinterLinInteraction*>(cBase);
                Mdouble unloadingStiffness = c->getUnloadingStiffness();
                c->setMaxOverlap(c->getMaxOverlap()
                                 *(unloadingStiffness-oldLoadingStiffness)
                                 /(unloadingStiffness-particleSpecies_->getLoadingStiffness())
                );
            }
        }
        oldTemperature = tempStep;
    }

    //--------------------------------------------------
    void actionsAfterTimeStep() override
    {
        setTemperature(getTemperature()); //With Temperature profile. i.e., dilatometer

        getNewVolume();

        //To compute the normal force
        for (auto i : interactionHandler)
        {
            scalarNormalForce += Vec3D::dot(i->getForce(),i->getNormal());
        }

        //To measure the mean coordination number.
        for (int i = particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        }
        meanCoordinationNumber /= particleHandler.getNumberOfObjects();

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
        Mdouble volSystem = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin()); //Todo!
        Mdouble volSintering = getNewVolume();//Todo!

        Mdouble volParticle = (4.0/3.0)*constants::pi*(std::pow(particleHandler.getMeanRadius(),3.0));
        Mdouble densityParticle = speciesHandler.getObject(0)->getDensity();

        Mdouble volTotalParticles = particleHandler.getVolume();
        Mdouble volumeFraction = volTotalParticles/volSystem;
        Mdouble massTotalParticles = particleHandler.getMass();

        Mdouble meanRadius = particleHandler.getMeanRadius();
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;

        Mdouble volParticlesPlusVoidInSystem = getNewVolume();//constants::pi*(std::pow(getXMax()/2,2.0))*getSurfaceHeight(); //Todo!
        Mdouble bulkDensity = massTotalParticles/(volParticlesPlusVoidInSystem);//Todo!

        Mdouble X = getMeanRelativeContactRadius()/2;
        Mdouble rho_p = 1.0/(1 + std::pow(1 - X,3.0));

        os << getTime() //[1] Current time
           << " " << getTimeMax() //[2] Max time
           << " " << volSystem //[3] Volume of the system
           << " " << getZMax() //[4] Height of the system
           << " " << particleHandler.getNumberOfObjects() //[5] Number of particles inserted
           << " " << particleHandler.getMeanRadius() // [6] Mean radius of a particle
           << " " << volParticle //[7] Mean volume of a particle
           << " " << densityParticle //[8] Mean density of a particle
           << " " << volTotalParticles //[9] Volume of particles inserted
           << " " << volumeFraction //[10] Volume fraction
           << " " << massTotalParticles //[11] total mass of particle inserted
           << " " << getNewVolume() // [12] Max height reached by particles.
           << " " << getKineticEnergy() // [13] kinetic energy
           << " " << getElasticEnergy() //[14] Elastic energy
           << " " << scalarNormalForce //[15] Normal force
           << " " << meanCoordinationNumber // [16] meanCoordinationNumber
           << " " << getTemperature()//[17]Current temperature
           << " " << particleSpecies_->getLoadingStiffness()// [18] Loading stiffness
           << " " << particleSpeciesAtMeltingPoint_->getLoadingStiffness()// [19] Loading stiffness
           << " " << getMeanRelativeContactRadius()//[20] Neckgrowth
           << " " << meanOverlap//[21]
           << " " << meanContactRadius //[22]
           << " " << getThermalEnergy() //[23]
           << " " << volSintering //[24]
           << " " << bulkDensity //[25]
           << " " << rho_p // [26]
           << std::endl;
    }

    void printTime() const override
    {
        Mdouble volSystem = (getXMax()-getXMin())* (getYMax()-getYMin())*(getZMax()-getZMin());
        std::cout << "t= " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ",tmax= " <<  std::setprecision(3) << std::left << std::setw(3) << getTimeMax()
                  << ",Temperature= " << std::setprecision(3) << std::left << std::setw(6) << getTemperature()
                  << ",pRadius= " << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getMeanRadius()
                  << ",K1= " << std::setprecision(3) << std::left << std::setw(6) << particleSpecies_->getLoadingStiffness()
                  << ",Heat Energy= " << std::setprecision(3) << std::left << std::setw(6) << getThermalEnergy()
                  << ",square(delta/R)= " << std::setprecision(3) << std::left << std::setw(6) << getMeanRelativeContactRadius()
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "Densification_PA12_Aged1_2";

    Mdouble meanRadius = 33.0e-6; //[mm]//
    Mdouble deltaR = 1.0e-3; //Relate to thermal expasion coefficient
    //define domain size
    Mdouble numLayers = 6;
    Mdouble domainLength = 6*(2*meanRadius); //[mm]
    Mdouble domainWidth = 6*(2*meanRadius); //[mm]
    Mdouble domainDepth = numLayers*(2*meanRadius); //[mm]

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 930.0; //[mg/mm^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 1.602e9; //Young's modulus
    const Mdouble eta = 0.4; //Poisson ratio
    const Mdouble pDepth = 1.45; //Penetration depth
    const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);

    const Mdouble beta = 0.0001; //Todo!
    const Mdouble k1 = beta*E_reduced*r_effective;
    const Mdouble k2 = 5.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.01; //

    //----------------------------------------------
    //Virgin. deltaC = 2.8e-7, C1 = 1.71, surface tension 0.4
    //1stAged: 2.9e-7, C1= 0.82, surface tension 0.4
    //2ndAgeed: deltaC = 3.00e-7, C1 = 0.497, surface tension 0.4
    Mdouble setDeltaC = 2.9e-07 ;
    Mdouble setC1 = 0.82;

    //Thermal variables:
    Mdouble setHCapacity = 1185.0; //[J/(kgK)]
    Mdouble setTConductivity = 0.231;  //[W/(mK)]

    Mdouble startingTemp_ = 165.0;
    Mdouble gradientTemp = 40.0;
    Mdouble meltTemp = 180.0;
    Mdouble maxTemp = 180.0;
    Mdouble timeCooling = 1.8;

    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //-------->[Set Plastic Contribution]
    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
    particleSpecies.setDensity(density);
    particleSpecies.setSinterAdhesion(0.0001*k1);

    //-------------------
    //Sinter parameters
    particleSpecies.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
    particleSpecies.setComplianceZero(1/(2*E_Modulus)); //Book: Adhesive Particles
    particleSpecies.setSurfTension(0.04); // [N/m] or [J/m^2] surface energy density or surface tension.

    //Parameters to calibrate:
    particleSpecies.setFluidity(setC1);
    particleSpecies.setSeparationDis(setDeltaC);

    particleSpecies.setThermalConductivity(setTConductivity);
    particleSpecies.setHeatCapacity(setHCapacity);

    //----------------------------------------------
    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    Densification pb (meanRadius,deltaR, numLayers,particleSpecies,startingTemp_,gradientTemp,
                      meltTemp, maxTemp, timeCooling);

    //----------------------------------------------
    pb.setGravity({0,0,0}); //mm/mu s^2 [Take care]
    pb.setMin({0.0,0.0,0.0});
    pb.setMax({domainLength,domainWidth,domainDepth});
    //----------------------------------------------

    Mdouble TimeStep = particleSpecies.getCollisionTime(mass)/50;

    //----------------------------------------------
    pb.setSaveCount(10000);
    pb.setTimeStep(TimeStep);
//    pb.setTimeStep(0.000001);

    pb.setTimeMax(2.0);

    logger(INFO,"Time step: %", pb.getTimeStep());
    //----------------------------------------------
//    pb.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    pb.setParticlesWriteVTK(false);

    pb.setName(setFilename);
    pb.setXBallsAdditionalArguments("-solidf -v0");
    pb.setFileType(FileType::ONE_FILE);

    pb.removeOldFiles();
    pb.setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");

    pb.solve();

    std::cout << "Execute gnuplot 'load 'Densification_PA12_virgin.gnu' ' to view output" << std::endl;
    helpers::writeToFile(setFilename + ".gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "plot 'Densification_PA12_virgin.fstat' u 1:20 title 'DEM simulation'  with lines linestyle 2\n"
                         //                         "replot 'PA12_Aged1_3.370e-5.txt' u ($1):($2) w p ls 1 title 'Experimental data 3.37 um' \n"
                         //                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
                         //                         "replot 'PA12_2.584e-5.txt' u ($1):($2) w p ls 2 title 'Experimental data 2.584' \n"
//                         "replot 'PA12_2.928e-5.txt' u ($1):($2) w p ls 3 title 'Experimental data 2.982 ' \n"
                         "replot 'PA12_3.044e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.044 '\n "
                         "replot 'PA12_3.210e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 ' \n"
                         "replot 'PA12_3.216e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 '\n "
//                         "replot 'PA12_2.93e-5_Zhao.txt' u ($1):($2) w p ls 4 title 'Experimental data 2.93 Zhao '\n "

    );

    std::cout << "Execute gnuplot 'load 'Densification_PA12_virgin.gnu' ' to view output" << std::endl;
    helpers::writeToFile(setFilename + "2.gnu",
                         "set xlabel 'a(t)/R'\n"
                         "set ylabel 'Q [J]'\n"
                         "plot 'Densification_PA12_virgin.fstat' u 20:23 title 'Heat energy'  with lines linestyle 2"
    );

    return 0;

}