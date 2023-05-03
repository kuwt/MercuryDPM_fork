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
#include <Boundaries/PeriodicBoundary.h>


using constants::pi;
using mathsFunc::cubic;

class DiscreteHeatTransfer : public Mercury3D
{
private:

    Mdouble radius_ = 0.0;
    Mdouble deltaR_ = 0.0; //[check] Thermal expansion coefficient
    Mdouble numLayers_ = 1.0; //1.0; By default at least one layer.

    Mdouble startingTemp_ = 0.0;
    Mdouble temperatureSurface_ = 0.0;
    Mdouble temperatureBottom_ = 0.0;

    Vec3D setHeatPosition_;

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint_;

public:

    DiscreteHeatTransfer (Mdouble radius, Mdouble deltaR, Mdouble numLayers, ParticleSpecies& particleSpecies, Mdouble  startingTemp, Mdouble temperatureSurface,
            Mdouble temperatureBottom)
    {
        radius_ = radius;
        deltaR_ = deltaR;
        numLayers_ = numLayers;

        startingTemp_ =  startingTemp;
        temperatureSurface_ = temperatureSurface;
        temperatureBottom_ = temperatureBottom;

        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions
        //-------------------------------------------------
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

        const unsigned numColumns = 10;
        const unsigned numrows = 10;

        Vec3D initial;
        initial.setZero();

        std::array<Vec3D,numColumns*numrows> generator = {initial};

        int counter = 0;
        for(int j = 1; j<=numColumns; j++){

            for(int i = 1; i<=numColumns; i++){

                Vec3D toFill(2*i+1,-2*j-1,0);
                generator[counter] = {toFill};
                counter += 1;
            }
        }



//        std::array<Vec3D, 4> layerPos = {Vec3D(1, 1, 0), Vec3D(1, -1, 0), Vec3D(-1, -1, 0), Vec3D(-1, 1, 0)};
//        std::array<Vec3D, 9> layerPos = {Vec3D(1, 1, 0), Vec3D(1, -1, 0), Vec3D(-1, -1, 0), Vec3D(-1, 1, 0),
//                                         Vec3D(3, 1, 0), Vec3D(3, -1, 0), Vec3D(3, -3, 0), Vec3D(1, -3, -0), Vec3D(-1, -3, 0)};
//        std::array<Vec3D, 36> layerPos = {Vec3D(1, 1, 0), Vec3D(1, -1, 0), Vec3D(1, -3, 0), Vec3D(-1, -3, 0),
//                                          Vec3D(1, -5, 0), Vec3D(-1, -5, 0),
//                                          Vec3D(-1, -1, 0), Vec3D(-1, 1, 0),
//                                          Vec3D(3, 1, 0), Vec3D(3, -1, 0),
//                                          Vec3D(3, -3, 0), Vec3D(3, -5, 0),
//                                          Vec3D(5, 1, 0), Vec3D(5, -1, -0),
//                                          Vec3D(5, -3, 0), Vec3D(5, -5, 0),
//                                          Vec3D(7, 1, 0), Vec3D(7, -1, 0),
//                                          Vec3D(7, -3, 0), Vec3D(7, -5, 0),
//                                          Vec3D(7, -7, 0), Vec3D(5, -7, 0),
//                                          Vec3D(3, -7, 0), Vec3D(1, -7, 0),
//                                          Vec3D(-1, -7, 0),Vec3D(9, 1, 0),
//                                          Vec3D(9, -1, 0),Vec3D(9, -3, 0),
//                                          Vec3D(9, -5, 0),Vec3D(9, -7, 0),
//                                          Vec3D(9, -9, 0),
//                                          Vec3D(7, -9, 0),Vec3D(5, -9, 0),
//                                          Vec3D(3, -9, 0),Vec3D(1, -9, 0),
//                                          Vec3D(-1, -9, 0)};
//        double posFirst = 11*radius_;
        double posFirst = radius_;
//        double posFirst = 1.0;
        for (unsigned i = 0; i < numLayers_; i++) {
            Vec3D center = Vec3D(getXMin()-2*radius_,
                                 getYMin()-2*radius_,
                                 getZMin()+ i * 1.999 * radius_ + posFirst);
            for (auto pos : generator) {
                P.setPosition(center + radius_ * pos*1.0);
//                P.fixParticle();
                P.setTemperature(startingTemp_);
                P.setVelocity(Vec3D(0.0, 0.0, 0.0));
                particleHandler.copyAndAddObject(P);
//                / fix the last particles in space
                if (i == numLayers_ - 1) {
                    particleHandler.getLastObject()->fixParticle();
//                    /// give an initial kick to the last layer of particles
////                    particleHandler.getLastObject()->setPosition(P.getPosition() + (Vec3D(0.0, -0.05 * radius_, 0.0)));
                }
            }
        }

        std::cout << std::endl;
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
        std::cout << std::endl;

        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0));
        w.set({0,0,-1},{0,0,getZMin()});
        wallHandler.copyAndAddObject(w);
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
    void setTemperature(Mdouble temperatureSurface)
    {

        //add temperature above a certain height
        Mdouble lb_z = getSurfaceHeight() - 2.0 * particleHandler.getMeanRadius();
        Mdouble lb_x = setHeatPosition_.X;
        Mdouble lb_y = setHeatPosition_.Y;

        //heat particles in top layer - Heat Transfer
        for (BaseParticle *p : particleHandler) {

            if (p->getPosition().Z > lb_z)
            {
                auto *tp = dynamic_cast<ThermalParticle *>(p);

                tp->setTemperature(temperatureSurface); //Change the temperature of the particle
            }else if(p->getPosition().Z < 4.0 * particleHandler.getMeanRadius()){

                auto *tp = dynamic_cast<ThermalParticle *>(p);

                tp->setTemperature(temperatureBottom_); //Change the temperature of the particle
            }
        }

    }

    void positionHeatSource()
    {
        //This vector will change eventually.
        Vec3D getCenter;
        getCenter.X = (getXMax()-getXMin())/2.0;
        getCenter.Y = (getYMax()-getYMin())/2.0;
        getCenter.Z = (getZMax()-getZMin())/2.0;

        setHeatPosition_.X =  getCenter.X;
        setHeatPosition_.Y = getCenter.Y;
        setHeatPosition_.Z = getSurfaceHeight();
    }

    //--------------------------------------------------
    void actionsAfterTimeStep() override
    {
        positionHeatSource();
        setTemperature(temperatureSurface_); //With Temperature profile. i.e., dilatometer

    }

    //--------------------------------------------------
    double getInfo(const BaseParticle& p) const override
    {
        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-startingTemp_)/(temperatureSurface_-startingTemp_);
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
           << " " << temperatureSurface_ //2
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
           << " " << getTotalMomentum() //17
           << " " << rho_p //18
           << std::endl;
    }

    // Write requested output to the ene file
    void writeEneTimeStep(std::ostream &os) const override
    {
        if (eneFile.getCounter() == 1)
        {
            os << "time Position Temperature NumParticles\n";
        }

        //get positions at the center of domain
        Mdouble position_x = (getXMax()-getXMin())/2;
        Mdouble position_y = (getYMax()-getYMin())/2;
        Mdouble position_z = (getZMax()-getZMin());


        Mdouble particleRadius = particleHandler.getMeanRadius();

        for (BaseParticle *p : particleHandler) {

            if (p->getPosition().X > (position_x - particleRadius) and p->getPosition().X < (position_x + particleRadius))
            {
//                if (p->getPosition().Y > (position_y))
//                {
                auto *tp = dynamic_cast<ThermalParticle *>(p);
                os << getTime() << ' '
                   << p->getPosition().Z << ' '
                   << tp->getTemperature() << ' '
                   << particleHandler.getNumberOfObjects() << '\n';
//                }
            }
        }

    }

    void printTime() const override
    {
        Mdouble volSystem = (getXMax()-getXMin())* (getYMax()-getYMin())*(getZMax()-getZMin());
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax=" <<  std::setprecision(3) << std::left << std::setw(3) << getTimeMax()
                  << ", Temperature=" << std::setprecision(3) << std::left << std::setw(6) << temperatureSurface_
                  << ", Heat Energy=" << std::setprecision(3) << std::left << std::setw(6) << getThermalEnergy()
                  << ", Boundary Heigth" << std::setprecision(3) << std::left << std::setw(6) << getZMax()
                  << ", Boundary Width" << std::setprecision(3) << std::left << std::setw(6) << getXMax()

                //                  << ", square(delta/R) " << std::setprecision(3) << std::left << std::setw(6) << getMeanRelativeContactRadius()
                  << std::endl;
    }

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "HeatTransferGrid2";

    Mdouble meanRadius = 100.0e-6; //[mm]//
    Mdouble deltaR = 1.0e-4; //Relate to thermal expasion coefficient
    Mdouble numLayers = 50;
    unsigned numParticles = 5;
    //define domain size
    Mdouble domainLength = numParticles*(2*meanRadius); //[mm]
    Mdouble domainWidth = numParticles*(2*meanRadius); //[mm]
    Mdouble domainDepth = numLayers*(2*meanRadius); //[mm]

    //----------------------------------------------
    //define properties of pcd article-particle contacts
    Mdouble density = 0.001; //[mg/mm^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 1.650; //Young's modulus
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble pDepth = 0.001; //Penetration depth
    const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);

    const Mdouble k1 = 0.001*E_reduced*r_effective;
    const Mdouble k2 = 5.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    //----------------------------------------------
    //Thermal variables:
    //PS
    Mdouble setHCapacity = 1185.0; //[J/(kgK)]
    Mdouble setTConductivity = 1.0;  //[W/(mK)]

    Mdouble startingTemp_ = 2.5;

    Mdouble temperatureSurface = 50.0;
    Mdouble temperatureBottom = 10.0;

    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //-------->[Set Plastic Contribution]
    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
    particleSpecies.setDensity(density);

    //-------------------
    particleSpecies.setThermalConductivity(setTConductivity);
    particleSpecies.setHeatCapacity(setHCapacity);

    //----------------------------------------------
    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    DiscreteHeatTransfer pb (meanRadius,deltaR,numLayers,particleSpecies,startingTemp_,
            temperatureSurface, temperatureBottom);

    //----------------------------------------------
    pb.setGravity({0,0,0.0}); //mm/mu s^2 [Take care]
    pb.setMin({0.0,0.0,0.0});
    pb.setMax({domainLength,domainWidth,domainDepth});
    //----------------------------------------------

    Mdouble TimeStep = particleSpecies.getCollisionTime(mass)/50;

    //----------------------------------------------
    pb.setSaveCount(10);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(0.011);

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
                         "plot 'HeatTransferGrid.fstat' u 1:12 title 'DEM simulation'  with lines linestyle 2\n"
    );

    helpers::writeToFile(setFilename + "2.gnu",
                         "set xlabel 'a(t)/R'\n"
                         "set ylabel 'Q [J]'\n"
                         "plot 'HeatTransferGrid.fstat' u 12:13 title 'Heat energy'  with lines linestyle 2"
    );

    return 0;
}