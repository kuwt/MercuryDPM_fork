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

#include "Math/Helpers.h"
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "Coupling/SurfaceCoupledSolidProblem.h"

/**
 * Define a coupled problem
 */
class CoupledBeam : public SurfaceCoupledSolidProblem<RefineableQDPVDElement<3, 2>>
{
public:

    CoupledBeam () {
        //set name
        setName("CoupledBeamRolling");
        //remove existing output files
        removeOldFiles();

        // setup steps
        setupOomph();
        setupMercury();

        // equalise time step
        if (getOomphTimeStep()>getTimeStep()) {
            setOomphTimeStep(getTimeStep());
        } else {
            setTimeStep(getOomphTimeStep());
        }

        // setup time
        setTimeMax(2.22);
        setSaveCount(4);
        logger(INFO,"timeMax: %, nTimeSteps %", getTimeMax(), getTimeMax()/getTimeStep());

        linear_solver_pt()->disable_doc_time();
        //disable_info_in_newton_solve();

        // Solve the  problem
        writeToVTK();
        solveSurfaceCoupling();
        saveSolidMesh();
    }

    void setupOomph() {
        setElasticModulus(1e6);
        setDensity(2500);
        setSolidCubicMesh(7, 5, 5, 0, 0.16, 0, 0.04, 0, 0.01);
        //pin at xmin and xmax
        pinBoundaries({Boundary::X_MIN, Boundary::X_MAX});
        // set time step
        //set time scale of oscillation
        //https://vlab.amrita.edu/?sub=3&brch=175&sim=1080&cnt=1
        double b = 0.04, d = 0.01, l = 0.08;
        double inertia = b*d*d*d/12;
        double omega = 1.875*1.875*sqrt(getElasticModulus()*inertia/(getDensity()*b*d*std::pow(l,4)));
        double timeScale = 2*constants::pi/omega; //0.06
        setOomphTimeStep(0.04*timeScale);
        logger(INFO,"Oomph time step %", getOomphTimeStep());
        //set dissipation
        //double criticalDissipation = 0.1 * sqrt(elasticModulus_ * density_ / lengthScale); //TW
        //setDissipation(criticalDissipation);
        //logger(INFO,"Adding dissipation %",criticalDissipation);
        //setBodyForceAsGravity();
        setNewtonSolverTolerance(3e-8);
        prepareForSolve();
        solveSteady();
        getBeamDeflection();
        // couple boundaries
        coupleBoundary(Boundary::Z_MAX);
    }

    void setupMercury() {
        //set domain
        std::array<double, 3> min;
        std::array<double, 3> max;
        getDomainSize(min, max);
        setDomain(min, max);
        // add elastic species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(1000);
        const double radius = 5e-3;
        const double mass = species->getMassFromRadius(radius);
        // stiffness such that k*(0.1*r)=.5*m*v^2
        double velocity = 0.05;
        double stiffness = mass*velocity*velocity/(0.01*radius*radius);
        species->setStiffness(stiffness);
        logger(INFO,"Species: %", *species);
        // add particle
        SphericalParticle p(species);
        p.setRadius(radius);
        double overlap = mass*9.8/species->getStiffness();
        p.setPosition({getXMin()+radius,getYCenter(),getZMax()+p.getRadius()-overlap});
        //p.setPosition({getXCenter(),getYMax()+p.getRadius(),getZCenter()});
        p.setVelocity({velocity,0,0});
        auto particle = particleHandler.copyAndAddObject(p);
        logger(INFO,"Particle: %", *particle);
        // set time step
        double tc = species->getCollisionTime(mass);
        setTimeStep(0.05*tc);
        logger(INFO,"Mercury time step %", getTimeStep());
        // set output file properties
        setParticlesWriteVTK(true);
        setWallsWriteVTK(true);
        setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        restartFile.writeFirstAndLastTimeStep();
        // add gravity
        setGravity({0,0,-9.8});
    }

    /// output file stream
    std::ofstream out;

    /// Write header of output file
    void actionsBeforeSolve() {
        helpers::writeToFile(
                getName()+".gnu",
                "set xlabel 'time'\n"
                "set ylabel 'kinetic and potential energy'\n"
                "p 'CoupledBeamRolling.out' u 1:($3+$4+$5) w lp t 'oomph-lib', '' u 1:(($6+$7+$8)) w lp t 'MercuryDPM', '' u 1:($3+$4+$5+$6+$7+$8) w lp t 'total'");
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
        xi[0] = 0.5*(max[0]+min[0]);
        xi[1] = 0.5*(max[1]+min[1]);
        xi[2] = 0.5*(max[2]+min[2]);
        return getDeflection(xi,2);
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
