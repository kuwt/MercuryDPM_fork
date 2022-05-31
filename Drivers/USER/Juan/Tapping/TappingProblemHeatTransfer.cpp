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

#include <Mercury3D.h>
#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>
//#include <Species/LinearViscoelasticSpecies.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Boundaries/DeletionBoundary.h>
using constants::pi;

class TapHeatTransfer: public Mercury3D
{
private:
    Mdouble cylinderRadius;
    Mdouble cylinderHeight;
    Mdouble radius;
    Mdouble bulkVolume;
    unsigned NumberOfTaps;

//    LinearViscoelasticSpecies* particleSpecies; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies;

public:
    //Constructor
    TapHeatTransfer(unsigned NumberOfTaps_, Mdouble radius_, ParticleSpecies& particleSpecies_)
    {
        NumberOfTaps = NumberOfTaps_;
        radius = radius_;
        logger(INFO,"Name: %",getName());

        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies_); //particle-particle interactions
        //-------------------------------------------------
        particleSpecies = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));

        double mass = particleSpecies->getMassFromRadius(radius);
        setTimeStep(0.02*particleSpecies->getCollisionTime(mass));

        cylinderRadius = 6.0*radius;
        cylinderHeight = 4.0*cylinderRadius;

        setMax({cylinderRadius,cylinderRadius,cylinderHeight});
        setMin({-cylinderRadius,-cylinderRadius,0});

        // apply a certain tapping velocity, upwards-fraction, and period
        const double a = radius; //tapping amplitude
        const double tr = 200*50*getTimeStep(); //time to rise (and relax)
        const double tf = sqrt(2.0*a/9.8); //time to fall
        const double tp = tr+tf; //tapping period

        logger(INFO,"Tapping period % ms, rise fraction %, tapping amplitude % mm", tp*1000, tr/tp, a*1000);
        auto tapPosition = [a, tr, tp] (double time) {
            const Mdouble t = fmod(time,tp)/tr; // scaled time, rising from 0 to 1 while the plate rises
            const Mdouble tFall = (t-1)*tr;
            return Vec3D(0,0,t<=1 ? a*t : std::max(0.0, a-0.5*9.8*tFall*tFall));
        };

        setTimeMax(NumberOfTaps*tp);
        restartFile.setSaveCount(tp/getTimeStep());

        //Walls
        InfiniteWall base;
        base.setSpecies(speciesHandler.getObject(0));
        base.set({0,0,-1},{0,0,0});
        base.setPrescribedPosition(tapPosition);
        wallHandler.copyAndAddObject(base);


        setName("TappingProblemHeatTransfer");
        logger(INFO,"Name: %",getName());

        AxisymmetricIntersectionOfWalls casing;
        casing.setSpecies(speciesHandler.getObject(0));
        casing.setAxis({0,0,1});
        casing.addObject({1,0,0},{cylinderRadius,0,0});
        casing.addObject({0,0,1},{0,0,0});
        casing.setPrescribedPosition(tapPosition);
        wallHandler.copyAndAddObject(casing);

        //add particles
        SphericalParticle particle;
        particle.setSpecies(speciesHandler.getObject(0));
        particle.setRadius(radius);
        CubeInsertionBoundary insertion;
        insertion.setHandler(&boundaryHandler);

        bulkVolume = 0.58*getTotalVolume()*constants::pi/4.0;
//        insertion.setPSD(psd);
        insertion.setInitialVolume(bulkVolume);
        Vec3D max = getMax()/0.3;

        do {
            insertion.set(&particle, 100, getMin(), max, Vec3D(0, 0, 0), Vec3D(0, 0, 0));
            insertion.checkBoundaryBeforeTimeStep(this);
            max.Z *= 1.05;
            logger(INFO,"Inserted % particles, fraction %",particleHandler.getSize(),particleHandler.getVolume()/bulkVolume);
        } while (insertion.getVolumeOfParticlesInserted()<bulkVolume);

    }

    // Write requested output to the ene file
    void writeEneTimeStep(std::ostream &os) const override
    {
        Mdouble MassParticlesInBed = particleHandler.getMass();

        if (eneFile.getCounter() == 1)
        {
            os << "time com wall eneRatio volumeFraction\n";
        }
        double comZ = getCentreOfMass().Z;
        double wallZ = wallHandler.getObject(0)->getPosition().Z;
        double heightMin = particleHandler.getVolume()/(cylinderRadius*cylinderRadius*constants::pi);
        double volumeFraction = 0.5 * heightMin / (comZ-wallZ);

        double volumeInserted = volumeFraction*getTotalVolume();
        double massInserted = particleSpecies->getDensity()*volumeInserted;

        double densityInserted = massInserted/volumeInserted;

        os  << getTime() << ' ' //[1]
            << comZ << ' ' //[2]
            << wallZ << ' ' //[3]
            << getKineticEnergy()/getElasticEnergy() << ' '//[4]
            << volumeFraction << ' ' //[5]
            << volumeInserted << ' ' //[6]
            << getTotalVolume() << ' ' //[7]
            << massInserted << ' ' //[8]
            << densityInserted << ' ' //[9]
            << particleHandler.getNumberOfObjects() << '\n'; //[10]
    }

    // Also write the ene information to the screen

    void printTime() const override {
        double comZ = getCentreOfMass().Z;
        double wallZ = wallHandler.getObject(0)->getPosition().Z;
        double heightMin = particleHandler.getVolume()/(cylinderRadius*cylinderRadius*constants::pi);
        double volumeFraction = 0.5 * heightMin / (comZ-wallZ);
        logger(INFO, "time % com % wall % ene % volFrac % %", getTime(), comZ, wallZ, getKineticEnergy()/getElasticEnergy(), volumeFraction, 0.5 * heightMin);
        //writeEneTimeStep(std::cout);
    }

};

//Main Function
int main(int argc, char* argv[])
{
    //Set problem parameters:
    unsigned NumberOfTaps = 20;
    Mdouble radius = 0.05;

    //create a contact law for particle-particle interactions
//    LinearViscoelasticSpecies particleSpecies;
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;


    particleSpecies.setDensity(1000);
    particleSpecies.setPlasticParameters(2e5, 10.0*2e5, 0.0, 0.01);
//    particleSpecies.setStiffness(2e5);
    particleSpecies.setDissipation(50);

    //----------------------------------------------
    //Thermal variables:
    //PS
    Mdouble setHCapacity = 1185.0; //[J/(kgK)]
    Mdouble setTConductivity = 0.35;  //[W/(mK)]

    Mdouble startingTemp_ = 25.0;

    Mdouble temperatureSurface = 500.0;
    Mdouble temperatureBottom = 100.0;

//        s.setThermalConductivity(setTConductivity);
//        s.setHeatCapacity(setHCapacity);


//        s.setSinterAdhesion(0.0001*2e5);
//        s.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
    Mdouble setDeltaC = 2.8e-07 ;
    Mdouble setC1 = 1.7;
    Mdouble compliance0 = 1.0/(2.0*2e5);
    Mdouble surfaceTension = 0.04;
//        s.setComplianceZero(compliance0); //Book: Adhesive Particles
//        s.setSurfTension(surfaceTension); // [N/m] or [J/m^2] surface energy density or surface tension.

    //Parameters to calibrate:
//        s.setFluidity(setC1);
//        s.setSeparationDis(setDeltaC);
    //----------------------------------------------

    //----------------------------------------------
    //Create a solver and run the commands in the constructor
    TapHeatTransfer problem(NumberOfTaps,radius,particleSpecies);

    problem.setGravity(Vec3D(0.0, 0.0, -9.81));
    problem.setXBallsAdditionalArguments("-v0 -solidf");

    helpers::writeToFile(problem.getName() + ".gnu",
                         "set xlabel 'Time [s]'\n"
                         "set ylabel 'COM_z [m]'\n"
                         "p [0.1:] '" + problem.getName() + ".ene' u 1:2 w l t ''");

    problem.setWallsWriteVTK(false);
    problem.setParticlesWriteVTK(false);

    problem.solve();
    return 0;
}
