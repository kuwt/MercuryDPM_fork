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

///This code tests our plastic force model, as published in Luding 2008.
class powdersAndgrains : public DPMBase {
public:

    explicit powdersAndgrains(Mdouble radius) {
        //-----------------
        //Global parameters
        std::string r = helpers::toString(radius);
        setName("PowdersAndGrainsPolyamide12_" + r);

        setFileType(FileType::ONE_FILE);
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
        const Mdouble pDepth = 1.35;

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
//        sf.setSurfTension(0.047); // [N/m] or [J/m^2] surface energy density or surface tension. Paper:Particle size limits for sintering polymer colloids. S Mazur.

        // Control parameters:
        sf.setSinterRate(0.115);
        sf.setFluidity(7.1e-6);
//        sf.setSeparationDis(9.8e-06);
        sf.setSeparationDis(0.0245*8*radius);

        sf.setSurfTension(0.047); // [N/m] or [J/m^2] surface energy density or surface tension.
        //http://www.surface-tension.de/solid-surface-energy.htm
        auto species = speciesHandler.copyAndAddObject(sf);

        setTimeMax(17);
        setTimeStep(0.00008 * collisionTime);
        setSaveCount(10000);

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

    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
//        logger(INFO," a/r %",sqrt(meanOverlap/meanRadius));
        return sqrt(meanOverlap/meanRadius);

    }

    void printTime() const override {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();

        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << ", r=" << std::setprecision(3) << std::left << std::setw(6)
                  << " Ene " << std::setw(12) << getKineticEnergy()/getElasticEnergy()
                  //                  << " X/a " << std::setw(12) << getMeanRelativeContactRadius()
                  << " square(delta/R) " << std::setw(12) << sqrt(meanOverlap/meanRadius)
                  << std::endl;
        //std::cout.flush();
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    powdersAndgrains sp0(3.39e-5);
    sp0.solve();

    //This helper is to see the Fn vs Overlap, rate of overlap.
    std::cout << "Execute gnuplot 'load 'PowdersAndGrains.gnu' ' to view output" << std::endl;
    helpers::writeToFile("PowdersAndGrains.gnu",
                         "set xlabel 'displacement'\n"
                         "set ylabel 'force'\n"
                         "plot 'PowdersAndGrainsPolyamide12_3.39e-05.fstat' u 1:9 w lp\n"
    );

    // time, i, j, x, y, z, delta, deltat, fn, ft, nx, ny, nz, tx, ty, tz

    //This helper is to see the neck growth
    std::cout << "Execute gnuplot 'load 'PowdersAndGrains2.gnu' ' to view output" << std::endl;
    helpers::writeToFile("PowdersAndGrains2.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "plot 'PowdersAndGrainsPolyamide12_3.39e-05.fstat' u ($1):(sqrt($7/3.39e-5)) title 'DEM simulation'  with lines linestyle 2\n"
                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
                         "replot 'PA12_Aged1_3.39e-5.txt' ($1):($2) w p ls 1 title 'Experimental data'"

    );

    return 0;

}
