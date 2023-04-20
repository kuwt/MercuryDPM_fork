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

// [!Headers]
#include "Mercury3D.h"
#include "Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h"
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/InfiniteWall.h>
#include "Walls/AxisymmetricIntersectionOfWalls.h"

// [!Headers]

// [!Main Class]
class Sintering : public Mercury3D
{
public:

    explicit Sintering (Mdouble meanRadius,Mdouble initialTemperature, Mdouble finalTemperature)
            : initialTemperature_(initialTemperature), finalTemperature_(finalTemperature)
    {
        // Global parameters:
        setName("DensSinteringHeatTrasnfer");
        setGravity({0,0,-9.8});
        setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");
        setFileType(FileType::ONE_FILE);
        setParticlesWriteVTK(true);
        wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
        TemperatureFile_.open("Temperature.data");

        //define domain size
        Mdouble domainLength = 260e-6;
        Mdouble domainWidth = 260e-6;
        Mdouble domainDepth = 1100e-6;

        setMin({-0.5*domainLength,-0.5*domainWidth,0});
        setMax({0.5*domainLength,0.5*domainWidth,domainDepth});

        //define properties of particle-particle contacts
        Mdouble density = 930.0;
        Mdouble restitution = 0.1;

        particleSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        const Mdouble mass = particleSpecies->getMassFromRadius(meanRadius);
        particleSpecies->setDensity(density);

        // ------------------------
        // Thermal properties:
        // https://designerdata.nl/materials/plastics/thermo-plastics/polyamide-12
        //c specific heat capacity
        particleSpecies->setHeatCapacity(1185.0); //J/kg/K
        //K thermal conduction
        particleSpecies->setThermalConductivity(0.231); //W/m/K
        // ------------------------

        const Mdouble k1 = 0.007;
        const Mdouble pDepth = 0.15;
//        const Mdouble k1 = 0.01 * meanRadius; //[Stiffness depends on particle radius]
//        const Mdouble pDepth = 1.35;

        particleSpecies->setStiffnessAndRestitutionCoefficient(k1, restitution, mass);
//        particleSpecies->setPlasticParameters(k1, 10.0*k1, k1, pDepth);
        particleSpecies->setPlasticParameters(k1, 10.0*k1, k1*0.01, pDepth);

        // Sintering partameters
        particleSpecies->setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);

        particleSpecies->setSinterAdhesion(0.0016 * k1);

        const Mdouble YoungM = 1600.0e6; //[Pa] Young's Modulus for polyamide
        particleSpecies->setComplianceZero(1.0/(2.0*YoungM)); //Book: Adhesive Particles
//        particleSpecies->setSurfTension(0.047); //
        particleSpecies->setSurfTension(0.0047); //

        particleSpecies->setSinterRate(0.12);
        particleSpecies->setFluidity(9.3e-05);
        particleSpecies->setSeparationDis(8.1e-06);

        //Wall species
        wallSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        wallSpecies-> setDensity(2500.0e6);
//        wallSpecies->setLoadingStiffness(particleSpecies->getLoadingStiffness());
        wallSpecies->setPlasticParameters(1000*k1, 10000*k1, 0.0, 0.001);

        speciesHandler.getMixedObject(particleSpecies,wallSpecies)->mixAll(particleSpecies,wallSpecies);
    }

    void setupInitialConditions() override
    {
        //Walls:
        InfiniteWall baseWall;
        baseWall.setSpecies(wallSpecies);
        baseWall.set(Vec3D(0.0,0.0,-1.0), Vec3D(0.0,0.0,getZMin()));
        wallHandler.copyAndAddObject(baseWall);

        AxisymmetricIntersectionOfWalls sideWall;
        sideWall.setSpecies(wallSpecies);
        sideWall.setPosition(Vec3D((getXMin()+getXMax())/2.0,(getYMin()+getYMax())/2.0,0));
        sideWall.setOrientation(Vec3D(0, 0, 1));
        sideWall.addObject(Vec3D(1,0,0),Vec3D((getXMax()-getXMin())/2.0,0,0));
        wallHandler.copyAndAddObject(sideWall);
    }

    void createColumn (Mdouble radius, Mdouble temperatureBottom, Mdouble temperatureRest)
    {
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(radius);
        p.setPosition({0,0,radius});
        p.setTemperature(temperatureBottom);
        particleHandler.copyAndAddObject(p);

        Mdouble domainRadius = getXMax()/2;
        Mdouble domainHeight = getZMax();
        Mdouble volumeOfParticles = 0.0;

        Mdouble addVolume  = 0.5*(getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());
        Mdouble fillHeight=getZMin();

        while (addVolume>0) {
            Mdouble x = random.getRandomNumber(-domainRadius+p.getRadius(), domainRadius-p.getRadius());
            Mdouble y = random.getRandomNumber(-domainRadius+p.getRadius(), domainRadius-p.getRadius());
            Mdouble z = random.getRandomNumber(getZMin(), fillHeight);

            p.setPosition({x, y, z});

            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                addVolume -= p.getVolume();
                p.setRadius(radius);
                p.setTemperature(temperatureRest);
                if (particleHandler.getNumberOfObjects()%100==0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects()%1000==0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects()%10000==0) std::cout << addVolume << '\n';
            } else {
                fillHeight += 0.01*radius; //increase fill height (slowly to insert particles as low as possible)
            }
        }
//        for (int i=1; i<numParticles;i++)
//        {
//            p.setPosition({0,0,radius+2*i*radius});
//            p.setTemperature(temperatureRest);
//            particleHandler.copyAndAddObject(p);
//        }
//        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());


        std::cout << std::endl;
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
        logger(INFO,"Total Volume particles",volumeOfParticles);
    }

    //remove particles according to (max.) outflow rate
    void actionsAfterTimeStep() override
    {
        //make decisions based on settled state
        static ThermalParticle* heatedParticle = nullptr;

        if (getTime()<0.05) return;

        particleSpecies->setPenetrationDepthMax(0.3);

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

            if (getKineticEnergy()/getElasticEnergy() < 5e-5 )
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

    // display time and eneRatio every time the files are printed
    void printTime() const override
    {
        int p1 = 0;
        int p2 = 1;
//        int p3 = 2;

        std::cout
                << "t " << std::setprecision(3) << std::left << std::setw(8)
                << getTime()
                << " EneRatio " << std::setprecision(3) << std::left << std::setw(8)
                << getKineticEnergy()/getElasticEnergy() <<std::endl
                << " T0 " << std::setprecision(3) << std::left << std::setw(6)
                << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(p1))->getTemperature()
                << "At Position:" <<particleHandler.getObject(p1)->getPosition()
                <<std::endl
                << " T1 " << std::setprecision(3) << std::left << std::setw(6)
                << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(p2))->getTemperature()
                << "At Position:" <<particleHandler.getObject(p2)->getPosition()
                << std::endl;
//         << " T2 " << std::setprecision(3) << std::left << std::setw(6)
//         << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(p3))->getTemperature()
//           << "At Position:" <<particleHandler.getObject(p3)->getPosition()
//         << std::endl;
    }

    double getInfo(const BaseParticle& p) const override
    {
        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-initialTemperature_)/(finalTemperature_-initialTemperature_);
        //return dynamic_cast<const ThermalParticle &>(p).getTemperature();
    }

private:

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* wallSpecies; //pointer to the wall properties

    Mdouble initialTemperature_;
    Mdouble finalTemperature_;

    std::ofstream TemperatureFile_;
};


// [!Main function]
int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble meanRadius = 31.0e-6;
    //define Temperatures
    Mdouble initialTemperature = 165.0;//K
    Mdouble finalTemperature = 180.0; //K
    //define particles in column

    Sintering pb (meanRadius,initialTemperature, finalTemperature);

//    pb.createColumn(numberParticlesInColumn,meanRadius,finalTemperature,initialTemperature);
//    pb.createColumn(numberParticlesInColumn,meanRadius,finalTemperature,finalTemperature);
    pb.createColumn(meanRadius,initialTemperature,initialTemperature);

    Mdouble collisionTime = 0.0008;
    pb.setTimeStep(collisionTime/50);
    pb. setSaveCount(40);
    pb.setTimeMax(3.0);

    pb.removeOldFiles();
    pb.solve();

    helpers::writeToFile("DensSinteringHeatTrasnfer.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'x/a'\n"
                         "plot 'DensSinteringHeatTrasnfer.fstat' u ($1):(sqrt($7/36.0e-6))");

    return 0;
}
