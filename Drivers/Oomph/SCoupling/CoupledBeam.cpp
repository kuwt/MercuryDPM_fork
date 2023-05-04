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

#include "Math/Helpers.h"
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "Oomph/SCoupling/SCoupledSolidProblem.h"

/**
 * Define a coupled problem
 */
class CoupledBeam : public SCoupledSolidProblem<RefineableQDPVDElement<3, 2>>
{
public:
    
    CoupledBeam () {
        //set name
        setName("CoupledBeam");
        //remove existing output files
        removeOldFiles();
        
        // setup steps
        setupOomph();
        setupMercury();

        // setup time
        setTimeMax(200);
        setSaveCount(500);
        logger(INFO,"timeMax: %, nTimeSteps %", getTimeMax(), getTimeMax()/getTimeStep());

        linear_solver_pt()->disable_doc_time();
        //disable_info_in_newton_solve();

        // Solve the  problem
        writeToVTK();
        solveSurfaceCoupling();
        saveSolidMesh();
    }
    
    void setupOomph() {
        setElasticModulus(100e6);
        setDensity(1309);
        double h = 0.4;
        setSolidCubicMesh(30, 2, 2, 0, 30*h, 0, 2*h, 0, 2*h);
        //pin at xmin and xmax
        pinBoundary(Boundary::X_MIN);
        // solve
        setNewtonSolverTolerance(1e-4);
        prepareForSolve();
        solveSteady();
        getBeamDeflection();
        // couple boundaries
        coupleBoundary(Boundary::Z_MAX);
    }

    void setupMercury() {
        //set domain
        std::array<double, 3> min, max;
        getDomainSize(min, max);
        setDomain(min, max);
        // add elastic species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(2e4);
        const double radius = 0.1;
        species->setStiffness(getElasticModulus()*2*radius);
        logger(INFO,"Species: %", *species);
        // add particle
        SphericalParticle p(species);
        p.setRadius(radius);
        const double mass = species->getMassFromRadius(radius);
        double overlap = mass*9.81/species->getStiffness();
        p.setPosition({0.2,0.2,getZMax()+radius-overlap});
        p.setVelocity({0.1,0,0});
        auto particle = particleHandler.copyAndAddObject(p);
        logger(INFO,"Particle: %", *particle);
        p.setPosition({0.2,0.6,getZMax()+radius-overlap});
        particle = particleHandler.copyAndAddObject(p);
        logger(INFO,"Particle: %", *particle);
        // set time step
        double tc = species->getCollisionTime(species->getMassFromRadius(p.getRadius()));
        setTimeStep(0.05*tc);
        // set oomph and mercury time step equal
        setOomphTimeStep(getTimeStep());
        logger(INFO,"Mercury time step %", getTimeStep());
        logger(INFO,"Oomph time step %", getOomphTimeStep());
        // set output file properties
        setParticlesWriteVTK(true);
        wallHandler.setWriteVTK(true);
        setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        restartFile.writeFirstAndLastTimeStep();
        // add gravity
        setGravity({0,0,-9.81});
    }

    /// output file stream
    std::ofstream out;

    /// Write header of output file
    void actionsBeforeSolve() override {
        helpers::writeToFile(
                getName()+".gnu",
                "set xlabel 'time'\n"
                "set ylabel 'kinetic and potential energy'\n"
                "p 'CoupledBeam.out' u 1:($3+$4) w lp t 'oomph-lib', '' u 1:(($6+$7-6.54498e-07 )) w lp t 'MercuryDPM', '' u 1:($3+$4+$6+$7-6.54498e-07) w lp t 'total'");
        out.open(getName()+".out");
        out << "time deflection elasticEnergy kineticEnergy gravEnergy elasticEnergyInteractions kineticEnergyParticles gravEnergyParticles particlePosition\n";
    }

    /// Each time step, compute deflection, elastic, kinetic and gravitational energy, and write to output file
    void actionsBeforeOomphTimeStep() override {
        double mass, elasticEnergy, kineticEnergy;
        Vector<double> com(3), linearMomentum(3), angularMomentum(3);
        getMassMomentumEnergy(mass, com, linearMomentum, angularMomentum, elasticEnergy, kineticEnergy);
        static double gravEnergy0 = getOomphGravity()*mass*com[2];
        double gravEnergy = getOomphGravity()*mass*com[2]-gravEnergy0;
        static double gravEnergyParticle0 = getGravitationalEnergy();
        out << getOomphTime()
            << ' ' << getBeamDeflection()
            << ' ' << elasticEnergy
            << ' ' << kineticEnergy
            << ' ' << gravEnergy
            << ' ' << getElasticEnergy()
            << ' ' << getKineticEnergy()
            << ' ' << getGravitationalEnergy()-gravEnergyParticle0;
        for (const auto p: particleHandler) {
            out << ' ' << p->getPosition();
        }
        out << std::endl;
    }

    /**
     * Outputs deflection at midpoint
     */
    double getBeamDeflection() {
        std::array<double,3> min, max;
        getDomainSize(min, max);
        
        Vector<double> xi(3);
        xi[0] = max[0];
        xi[1] = 0.5*(max[1]+min[1]);
        xi[2] = 0.5*(max[2]+min[2]);
        double deflection = getDeflection(xi,2);
        logger(INFO, "Beam deflection at center (% % %) is %",
               xi[0], xi[1],xi[2], deflection);

        double mass, elasticEnergy, kineticEnergy;
        Vector<double> com(3), linearMomentum(3), angularMomentum(3);
        getMassMomentumEnergy(mass, com, linearMomentum, angularMomentum, elasticEnergy, kineticEnergy);
        logger(INFO, "mass %, linearMomentum %, angularMomentum %, elasticEnergy %, totalEnergy %",
               mass, linearMomentum, angularMomentum, elasticEnergy, elasticEnergy+kineticEnergy);

        return deflection;
    }
};

/**
 * Measure bag deformation loaded by body force
 */
int main()
{
    CoupledBeam problem;
    return 0;
}
