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


#include "DPMBase.h"
#include <iostream>
#include <vector>
#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Logger.h>

///This code tests the linear plastic viscoelastic [and friction] behavior of two particles in sintering.
class calibrationUnitTest : public DPMBase {
public:

    explicit calibrationUnitTest(Mdouble radius) {
        //-----------------nu
        //Global parameters
        std::string r = helpers::toString(radius);
        setName("LinearPlasticViscoElasticFrictionThermalTest" + r);

        setFileType(FileType::ONE_FILE);
        setSaveCount(10000);
        setTimeMax(7);
        //setParticlesWriteVTK(true);

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
        const Mdouble k1 = 0.01 * radius; //[Stiffness depends on particle radius]
        const Mdouble restitutionCoefficient = 0.1;
        const Mdouble pDepth = 0.95;

        const Mdouble YoungM = 2.95e9; //[Pa] Young's Modulus for polystyrene

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
//        sf.setSinterRate(0.16);
        sf.setSinterAdhesion(0.0016 * k1);

        //-------------------
        //Parameters to be used in contact growth driven by adhesive intersurface forces and accomodated by
        //viscoelastic deformation.
        sf.setComplianceZero(1/(2*YoungM)); //Book: Adhesive Particles
        sf.setSurfTension(0.047); // [N/m] or [J/m^2] surface energy density or surface tension. Paper:Particle size limits for sintering polymer colloids. S Mazur.
        sf.setFluidity(4.0e-6);
        sf.setSeparationDis(2.8e-5);

        auto species = speciesHandler.copyAndAddObject(sf);

        setTimeStep(0.00007 * collisionTime);
        //-------------------
        //Particle properties:
        ThermalParticle P0, P1;
        //SphericalParticle P0, P1;
        P0.setSpecies(species);
        P1.setSpecies(species);

        P0.setRadius(radius);
        P1.setRadius(radius);

        P0.setPosition(Vec3D(-(1 - 1e-15) * radius, 0, 0));
        P1.setPosition(-P0.getPosition());

        particleHandler.copyAndAddObject(P0);
        particleHandler.copyAndAddObject(P1);

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //Friction contribution
        const Mdouble friction =  0.7;
        Mdouble k_s = 0.2*sf.getLoadingStiffness();
        Mdouble k_r = 0.1*sf.getLoadingStiffness();
        Mdouble k_o = 0.1*sf.getLoadingStiffness();

        Mdouble miu_s = 1.0*friction;
        Mdouble miu_r = 0.1*friction;
        Mdouble miu_o = 0.1*friction;

        Mdouble gamma_s = 0.2*sf.getDissipation();
        Mdouble gamma_r = 0.05*sf.getDissipation();
        Mdouble gamma_o = 0.05*sf.getDissipation();

        sf.setSlidingStiffness(k_s);
        sf.setRollingStiffness(k_r);
        sf.setTorsionStiffness(k_o);

        sf.setSlidingFrictionCoefficient(miu_s);
        sf.setRollingFrictionCoefficient(miu_r);
        sf.setTorsionFrictionCoefficient(miu_o);

        sf.setSlidingDissipation(gamma_s);
        sf.setRollingDissipation(gamma_r);
        sf.setTorsionDissipation(gamma_o);

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    }

    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
        logger(INFO," a/r %",sqrt(2.0*meanOverlap/meanRadius));
        return sqrt(2.0*meanOverlap/meanRadius);

    }

    void printTime() const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << ", r=" << std::setprecision(3) << std::left << std::setw(6)
                  << " Ene " << std::setw(12) << getKineticEnergy()/getElasticEnergy()
                  << " X/a " << std::setw(12) << getMeanRelativeContactRadius()
                  //<< ", Ekin=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
                  << std::endl;
        //std::cout.flush();
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    calibrationUnitTest sp0(3.216e-5);
    sp0.solve();

    //This helper is to see the Fn vs Overlap, rate of overlap.
    std::cout << "Execute 'gnuplot LinearPlasticViscoElasticFrictionThermalTest.gnu' to view output" << std::endl;
    helpers::writeToFile("LinearPlasticViscoElasticFrictionThermalTest.gnu",
                         "set xlabel 'displacement'\n"
                         "set ylabel 'force'\n"
                         "plot 'LinearPlasticViscoElasticFrictionThermalTest3.216e-05.fstat' u 7:9 w lp\n"
    );

    //This helper is to see the neck growth
    std::cout << "Execute 'gnuplot LinearPlasticViscoElasticFrictionThermalTest.gnu' to view output" << std::endl;
    helpers::writeToFile("LinearPlasticViscoElasticFrictionTest.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "plot 'LinearPlasticViscoElasticFrictionThermalTest3.216e-05.fstat' u ($1):(sqrt($7/3.216e-5)) title 'DEM simulation'  with lines linestyle 2\n"
                         "replot 0.3*x**0.5 title 'Frenkel' with lines linestyle 1"
    );

    return 0;

}
