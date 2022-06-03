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
#include "Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h"
#include "Walls/InfiniteWall.h"

// [!Main Class]
class heatTransferGranular : public Mercury3D
{
public:

    explicit heatTransferGranular (int numberParticlesInColumn, Mdouble meanRadius,Mdouble initialTemperature, Mdouble finalTemperature,Mdouble setHCapacity, Mdouble setTConductivity)
    {
        pRadius_ = meanRadius;
        NumParticles_ = numberParticlesInColumn;

        initialTemperature_ = initialTemperature;
        finalTemperature_ = finalTemperature;

        // Global parameters:
        setName("HeatTransferGranularMedium");
        setGravity({0.0,0.0,-9.81});

//        setGravity({0,0,0.0});
        setParticlesWriteVTK(true);
        setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");
        TemperatureFile_.open("Temperature.data");

        //-------------------
        //Boundary parameters
        //define domain size
        Mdouble domainLength = 2.0*pRadius_;
        Mdouble domainWidth = 2.0*pRadius_;
        Mdouble domainDepth = NumParticles_*2.0*pRadius_;

        setMin({-0.5*domainLength,-0.5*domainWidth,0});
        setMax({0.5*domainLength,0.5*domainWidth,domainDepth});

        //-------------------
        //Setup parameters:
        // Polystyrene
        //define properties of particle-particle contacts
        Mdouble density = 1000.0;
        Mdouble restitution = 0.1;
        const Mdouble pDepth = 0.01;

        const Mdouble k1 = 0.01;

        particleSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        const Mdouble mass = particleSpecies->getMassFromRadius(pRadius_);
        particleSpecies->setDensity(density);

        particleSpecies->setRestitutionCoefficient(restitution,mass);
        particleSpecies->setPlasticParameters(k1, 10*k1, 0.0, pDepth);
        particleSpecies->setSinterAdhesion(0.001*k1);

        // ------------------------
        // Thermal properties:
        particleSpecies->setThermalConductivity(setTConductivity);
        particleSpecies->setHeatCapacity(setHCapacity);
        // ------------------------

        // ------------------------
        //Wall species
        wallSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        wallSpecies-> setDensity(2500.0e6);
//        wallSpecies->setLoadingStiffness(particleSpecies->getLoadingStiffness());
        wallSpecies->setPlasticParameters(1000*k1, 10000*k1, 0.0, 0.001);

        // ------------------------
        speciesHandler.getMixedObject(particleSpecies,wallSpecies)->mixAll(particleSpecies,wallSpecies);
    }

    void setupInitialConditions() override
    {
        logger(INFO,"Creating a base wall at z=0");
        InfiniteWall w;
        w.setSpecies(wallSpecies);
        w.set({0,0,-1},{0,0,getZMin()});
        wallHandler.copyAndAddObject(w);
    }

    void createColumn ()
    {
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(pRadius_);
        p.setPosition({0,0,pRadius_});
        p.setTemperature(initialTemperature_);
        particleHandler.copyAndAddObject(p);

        for (int i=1; i<NumParticles_;i++)
        {
            p.setPosition({0,0,pRadius_+2*i*pRadius_});
            p.setTemperature(initialTemperature_);
            particleHandler.copyAndAddObject(p);
        }
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());
    }

    //remove particles according to (max.) outflow rate
    void actionsAfterTimeStep() override
    {
        //make decisions based on settled state
        static ThermalParticle* heatedParticle = nullptr;

//        if (getTime()<1.0) return;

        TemperatureFile_ << (getTime()-1.0) <<"\t";

        for (int i=0;i<particleHandler.getNumberOfObjects();i++)
        {
            double particleTemperture =dynamic_cast<const ThermalParticle*>(particleHandler.getObject(i))->getTemperature();
            TemperatureFile_ << particleTemperture << "\t";

        }
        TemperatureFile_ << "\n";

        if (heatedParticle == nullptr)
        {
            // Functions after this statement only get executed every 100th time step (if counter==100)
            static unsigned counter = 0;
            if (++counter != 100) return;
            else counter = 0;

            if (getKineticEnergy() < 5e-10 )
            {
                //set heated particle
                Mdouble minDistance = std::numeric_limits<Mdouble>::infinity();
                for (auto p : particleHandler) {
                    Mdouble distance = Vec3D::getDistance(p->getPosition(),Vec3D(0,0,getZMax()));
                    if (distance < minDistance) {
                        minDistance = distance;
                        heatedParticle = dynamic_cast<ThermalParticle*>(p);
                    }
                }
                logger(INFO, "Starting to heat the material");
            }
        } else {
            heatedParticle->setTemperature(finalTemperature_);
        }
    }

    //System
    Mdouble getThermalEnergy() const //
    {
        Mdouble eneEnergy;

        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;

            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            eneEnergy += (tp->getTemperature()-initialTemperature_) * tp->getMass() * particleSpecies->getHeatCapacity();
//            eneEnergy += deltaTemperature * tp->getMass() * particleSpecies->getHeatCapacity();
        }
        return eneEnergy;
    }

    void writeFstatHeader(std::ostream &os) const override
    {
//        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
//        Mdouble meanRadius = particleHandler.getMeanRadius();
//        Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;
        Mdouble deltaT = (dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getTemperature() - initialTemperature_);
        Mdouble eneThermal_i = deltaT * dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getMass()* particleSpecies->getHeatCapacity();

        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = 2.0*particleHandler.getMeanRadius();
        Mdouble meanContactRadius = (sqrt(meanOverlap/meanRadius)*meanRadius);

//        Mdouble distancePropaga = (meanContactRadius/eneThermal_i)*particleSpecies->getThermalConductivity()*deltaT;

        os << getTime() //1
           << " " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(1))->getTemperature() //2
           << " " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getTemperature() //3
           << " " << getThermalEnergy()//4[J]
           << " " << getThermalEnergy()/(getTime()) //5 [W/s]
           << " " << eneThermal_i //6
           << " " << getKineticEnergy() //7
           << " " << getElasticEnergy() //8
           << " " << meanRadius
           << std::endl;
    }

    double getInfo(const BaseParticle& p) const override
    {
        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-initialTemperature_)/(finalTemperature_-initialTemperature_);
        //return dynamic_cast<const ThermalParticle &>(p).getTemperature();
    }

    // display time and eneRatio every time the files are printed
    void printTime() const override
    {
        int p1 = 0;
        int p2 = 1;
//        int p3 = 2;

        std::cout
                << "t " << std::setprecision(3) << std::left << std::setw(8)
                << getTime()
                << " EneRatio " << getThermalEnergy()
                << ", Temperature at top " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(1))->getTemperature()
                << ", Temperature at bottom " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getTemperature()
                << ", Thermal energy = " << getThermalEnergy()
                << std::endl;
    }


private:

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* wallSpecies; //pointer to the wall properties

    Mdouble initialTemperature_;
    Mdouble finalTemperature_;

    std::ofstream TemperatureFile_;

    Mdouble pRadius_ = 1.0;
    int NumParticles_ = 1;
};


// [!Main function]
int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble meanRadius = 6.0e-5;
    //define Temperatures
    Mdouble initialTemperature = 0.0;//C
    Mdouble finalTemperature = 240.0; //C

    Mdouble setHCapacity = 1320.0; //[J/(kgK)]
    Mdouble setTConductivity = 0.167;  //[W/(mK)]

    //define particles in column
   int numberParticlesInColumn = 10;

    heatTransferGranular pb (numberParticlesInColumn,meanRadius,initialTemperature, finalTemperature,setHCapacity,setTConductivity);

//    pb.createColumn(numberParticlesInColumn,meanRadius,finalTemperature,initialTemperature);
//    pb.createColumn(numberParticlesInColumn,meanRadius,finalTemperature,finalTemperature);
    pb.createColumn();

    Mdouble collisionTime = 0.002;
    pb.setTimeStep(collisionTime/50);
    pb. setSaveCount(40);
    pb.setTimeMax(3.0);

    pb.removeOldFiles();
    pb.solve();

    helpers::writeToFile("SinteringPlusHeatTransfer.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'x/a'\n"
                         "plot 'SinteringPlusHeatTransfer.fstat' u ($1):(sqrt($7/36.0e-6))");

    return 0;
}
