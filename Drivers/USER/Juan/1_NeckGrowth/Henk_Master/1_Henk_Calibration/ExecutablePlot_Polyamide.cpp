//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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
#include <iostream>
#include <vector>
#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Logger.h>

///This code tests the sintering of two particles following three different contact interactions.

class calibrationSinterForceUnitTest : public DPMBase {
public:

    explicit calibrationSinterForceUnitTest(Mdouble radius,Mdouble tsintering,Mdouble materialC1, Mdouble separationDc) {

        //-----------------
        //Global parameters
        std::string r = helpers::toString(radius);
        std::string sint = helpers::toString(tsintering);
        std::string mtc1 = helpers::toString(materialC1);
        std::string sDc = helpers::toString(separationDc);

        setName("Test_" + sint + "_" + mtc1 + "_" + sDc);

        //  setFileType(FileType::ONE_FILE);
        dataFile.setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::NO_FILE);
        fStatFile.setFileType(FileType::ONE_FILE);
        eneFile.setFileType(FileType::NO_FILE);
        //setParticlesWriteVTK(true);
        //wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
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
        const Mdouble density = 1000;
        const Mdouble k1 = 10.0 * radius; //[Stiffness depends on particle radius]
        const Mdouble restitutionCoefficient = 0.1;
        const Mdouble pDepth = 1.4;

        const Mdouble YoungM = 1.65e9; //[Pa] Young's Modulus for polyamide

        //-------------------
        //Species:
        ThermalSinterLinFrictionReversibleAdhesiveSpecies sf;
        //SinterSpecies sf;
        sf.setHandler(&speciesHandler);
        sf.setDensity(density);

        const Mdouble mass = sf.getMassFromRadius(radius);
        sf.setStiffnessAndRestitutionCoefficient(k1, restitutionCoefficient, mass);
        sf.setPlasticParameters(k1, 10.0 * k1, k1, pDepth);

        const Mdouble collisionTime = sf.getCollisionTime(mass);
        logger(INFO, "Collision time %", collisionTime);

        //-------------------
        //Sinter parameters
        sf.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
        sf.setSinterAdhesion(0.0016 * k1);

        //-------------------
        //Parameters to be used in contact growth driven by adhesive intersurface forces and accomodated by
        //viscoelastic deformation.
        sf.setComplianceZero(1/(2*YoungM)); //Book: Adhesive Particles
        sf.setSurfTension(0.041); // [N/m] or [J/m^2] surface energy density or surface tension. Paper:Particle size limits for sintering polymer colloids. S Mazur.

        //-------------------
        //Parameters that don't remain constant during the simulations
//        sf.setSinterRate(tsintering);
        sf.setFluidity(materialC1);
        sf.setSeparationDis(separationDc);

        auto species = speciesHandler.copyAndAddObject(sf);

        setTimeMax(20.0); //4

        setTimeStep(0.0006 * collisionTime);
        setSaveCount(1000);
        //setSaveCount(1000);
        //-------------------
        //Particle properties:
        ThermalParticle P0, P1;
        //SphericalParticle P0, P1;
        P0.setSpecies(species);
        P1.setSpecies(species);

        P0.setRadius(radius);
        P1.setRadius(radius);

        P0.setPosition(Vec3D(-(1 - 1.0e-5) * radius, 0, 0));
        P1.setPosition(-P0.getPosition());

        particleHandler.copyAndAddObject(P0);
        particleHandler.copyAndAddObject(P1);

    }

//    void writeFstatHeader(std::ostream &os) const override {
//
//        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
//        Mdouble meanRadius = particleHandler.getMeanRadius();
//
//        os
////            << getTime()
//            << " " << sqrt(meanOverlap/meanRadius)
//            << std::endl;
//    }

    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
        logger(INFO," a/r %",sqrt(meanOverlap/meanRadius));
        return sqrt(meanOverlap/meanRadius);

    }

//    void printTime() const override {
//        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
//                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
//                  << ", r=" << std::setprecision(3) << std::left << std::setw(6)
//                  << " Ene " << std::setw(12) << getKineticEnergy()/getElasticEnergy()
//                  << " X/a " << std::setw(12) << getMeanRelativeContactRadius()
//                  << std::endl;
//    }
    void printTime() const override {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();

        std::cout << " " << sqrt(meanOverlap/meanRadius)
                  << std::endl;
    }
};

int main(int argc, char *argv[])
{
    //The following parameters need to be inserted after generating the executable.
    //- Generate the executable: make Calibration_Polyamide
    //- Run the executable: ./Calibration_Polyamide 2.58e-5 0.170 5.40e-6 2.80e-5

    //Information:
    //argv[0] = Calibration_Polyamide
    //argv[1] = 2.58e-5 : Particle radius
    //argv[2] = 0.170 : Sintering time
    //argv[3] = 5.40e-6 : Material constant
    //argv[4] = 2.80e-5 : Separation distance

    //Particle radius
    Mdouble pRadius = atof(argv[1]);

    //Sintering time
    Mdouble tsintering = atof(argv[2]);

    //Material constant
    Mdouble materialC1 = atof(argv[3]);

    //Separation distance
    Mdouble separationDc = atof(argv[4]);

    calibrationSinterForceUnitTest sp0(pRadius,tsintering,materialC1,separationDc);

    sp0.solve();

    return 0;

}
