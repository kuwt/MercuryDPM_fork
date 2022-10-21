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

#include "DPMBase.h"
#include <iostream>
#include <vector>
#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Logger.h>

///This code tests our plastic force model, as published in Luding 2008.
class powdersAndgrains : public DPMBase {
public:

    explicit powdersAndgrains(Mdouble radius) {
        //-----------------
        //Global parameters
        std::string r = helpers::to_string(radius);
        setName("PowdersAndGrainsPolyamide12_" + r);

        setFileType(FileType::ONE_FILE);
        //setParticlesWriteVTK(true);

        //setWallsWriteVTK(FileType::MULTIPLE_FILES);

        setParticleDimensions(3);
        setSystemDimensions(3);
        setGravity(Vec3D(0.0, 0.0,0.0));

        //-------------------
        //Boundary parameters
        setXMax(2.0 * radius);
        setYMax(radius);
        setZMax(radius);
        setXMin(-getXMax());
        setYMin(-getYMax());
        setZMin(-getZMax());

        //-------------------
        //Setup parameters:
//        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltTemp- getTemperature())/varTemp));
        const Mdouble k1 = (0.05 * radius); //[Stiffness depends on particle radius]

        const Mdouble restitutionCoefficient = 0.1;
        const Mdouble pDepth = 1.45;
        const Mdouble YoungM = 2.95e9; //[Pa] Young's Modulus for polyamide

        //-------------------
        // Species:
        speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        particleSpecies = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));

        particleSpeciesAtMeltingPoint = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies);

        const Mdouble mass = particleSpecies->getMassFromRadius(radius);
        particleSpecies->setStiffnessAndRestitutionCoefficient(k1, restitutionCoefficient, mass);
        particleSpecies->setPlasticParameters(k1, 10.0 * k1, k1, pDepth);

        //-------------------
        // Sintering parameters
        particleSpecies->setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
        particleSpecies->setSinterAdhesion(0.0016 * k1);

        //-------------------
        particleSpecies->setComplianceZero(1/(2*YoungM));
        particleSpecies->setSurfTension(0.047);

        // Control parameters:
//        particleSpecies->setSinterRate(0.064);
        particleSpecies->setFluidity(3.7e-06);
        particleSpecies->setSeparationDis(1.6e-05);

        //-------------------
        //Particle properties:
        ThermalParticle P0, P1;
        //SphericalParticle P0, P1;
        P0.setSpecies(particleSpecies);
        P1.setSpecies(particleSpecies);

        P0.setRadius(radius);
        P1.setRadius(radius);

        P0.setPosition(Vec3D(-(1 - 1.0e-5) * radius, 0, 0));
        P1.setPosition(-P0.getPosition());

        particleHandler.copyAndAddObject(P0);
        particleHandler.copyAndAddObject(P1);

        // Integration time:
        const Mdouble collisionTime = particleSpecies->getCollisionTime(mass);
        setTimeMax(7);
        setTimeStep(0.00008 * collisionTime);
        setSaveCount(10000);
    }

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies;
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint;

    // Initial conditions
    void setupInitialConditions() override
    {
//        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltTemp- getTemperature())/varTemp));
        particleSpeciesAtMeltingPoint->setLoadingStiffness(particleSpecies->getLoadingStiffness());
    }

    // Before the program runs
    void actionsBeforeTimeStep() override
    {
        const Mdouble density = 1000;
        particleSpecies->setDensity(density);
    }

    // Defining temperature
    Mdouble getTemperature()const
    {
        Mdouble heatingTemperature = 150.0 + tempGradient*getTime();
//        Mdouble coolingTemperature = 20.0 + tempGradient*(getTimeMax() - getTime());
//        return std::min(std::min(heatingTemperature,coolingTemperature),maxTemp);
        return std::min(heatingTemperature,maxTemp);
    }


    void setTemperature(Mdouble temperature)
    {
        static Mdouble oldTemperature = 150.0;
        Mdouble deltaTemperature = temperature-oldTemperature;
        Mdouble factorRadius = 1.0-deltaR*deltaTemperature;

        particleSpecies->setDensity(mathsFunc::cubic(factorRadius)*particleSpecies->getDensity());

        for (BaseParticle* p : particleHandler)
        {
            if (p->getSpecies()==particleSpecies) {
                p->setRadius(factorRadius * p->getRadius());
                NewRadius = p->getRadius();
            }
        }

        //change species properties
        Mdouble stiffnessFactor = 0.5*(1+tanh((meltTemp-temperature)/varTemp));

        particleSpecies->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint->getLoadingStiffness());

        //change old temperature
        oldTemperature = temperature;
    }

    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
//        logger(INFO," a/r %",sqrt(meanOverlap/meanRadius));
        return sqrt(meanOverlap/meanRadius);

    }

    void actionsAfterTimeStep() override
    {
        setTemperature(getTemperature());
    }

    void printTime() const override {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();

//        for (const auto& q: particleHandler) {
//            auto p0 = dynamic_cast<ThermalParticle *>(q);
//            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
//            logger(INFO, "time % temperature % moltenLayer % radius % solidRadius % invMass % mass % density % volume %",
//                   getTime(), p0->getTemperature(), p0->getMoltenLayerThickness(),
//                   p0->getRadius(), p0->getSolidRadius(), p0->getInvMass(), p0->getMass(), p0->getSpecies()->getDensity(), p0->getVolume());
//        }

        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << ", Temp " << std::setw(6) << getTemperature()
                  << ", Density " << std::setw(3) << particleSpecies->getDensity()
                  << ", newRadius " << std::setw(3) << NewRadius
                  << ", loading stiffness " << std::setw(3) << particleSpecies->getLoadingStiffness()
//                  << ", New stiffness " << std::setw(3) <<NewStiffness
                  << ", square(delta/R) " << std::setw(12) << sqrt(meanOverlap/meanRadius)
                  << std::endl;
        //std::cout.flush();
    }

private:

    Mdouble maxTemp = 190.0; //[C]
    Mdouble meltTemp = 180.0; //[C]
    Mdouble varTemp = 30; //[C]
    Mdouble deltaR = 1.0e-4; //[check]
    Mdouble tempGradient = 25.0; //[C]
    Mdouble NewRadius=1.0;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble radius = 3.055e-5; // [micro meters]

    powdersAndgrains sp0(radius);
    sp0.removeOldFiles();
    sp0.solve();

    //This helper is to see the Fn vs Overlap, rate of overlap.
    std::cout << "Execute gnuplot 'load 'PowdersAndGrains.gnu' ' to view output" << std::endl;
    helpers::writeToFile("PowdersAndGrains.gnu",
                         "set xlabel 'displacement'\n"
                         "set ylabel 'force'\n"
                         "plot 'PowdersAndGrainsPolyamide12_3.055e-05.fstat' u 1:9 w lp\n"
    );

    // time, i, j, x, y, z, delta, deltat, fn, ft, nx, ny, nz, tx, ty, tz

    //This helper is to see the neck growth
    std::cout << "Execute gnuplot 'load 'PowdersAndGrains2.gnu' ' to view output" << std::endl;
    helpers::writeToFile("PowdersAndGrains2.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "plot 'PowdersAndGrainsPolyamide12_3.055e-05.fstat' u ($1):(sqrt($7/3.055e-5)) title 'DEM simulation'  with lines linestyle 2\n"
                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
                         "replot 'PA12_Aged1_3.055e-5.txt' ($1):($2) w p ls 1 title 'Experimental data'"

    );

    return 0;

}
