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
#include "SurfaceCoupledSolidProblem.h"

/**
 * Define a coupled problem
 */
class CoupledBeam : public SurfaceCoupledSolidProblem<RefineableQDPVDElement<3, 2>>
{
public:
    
    CoupledBeam () {
        //set name
        setName("CoupledBeamUnitTest");
        //remove existing output files
        removeOldFiles();
        
        // setup steps
        setupOomph();
        setupMercury();

        // setup time
        setTimeMax(4e-3);
        setSaveCount(4);
        logger(INFO,"timeMax: %, nTimeSteps %", getTimeMax(), getTimeMax()/getTimeStep());

        linear_solver_pt()->disable_doc_time();
        //disable_info_in_newton_solve();

        // Solve the  problem
        writeToVTK();
        solveSurfaceCoupling();
        saveSolidMesh();

        double velocity = particleHandler.getLastObject()->getVelocity().Z;
        //logger(INFO,"Final particle velocity %", velocity);
        helpers::check(velocity, 0.506165, 1e-6, "final particle velocity");
    }
    
    void setupOomph() {
        setElasticModulus(1e6);
        setDensity(2500);
        setSolidCubicMesh(7, 5, 5, 0, 0.16, 0, 0.04, 0, 0.01);
        //pin at xmin and xmax
        pinBoundaries({Boundary::X_MIN, Boundary::X_MAX});
        // set time step
        double lengthScale = 0.1;
        double timeScale = 0.38 * lengthScale * sqrt(density_ / elasticModulus_); //TW
        setOomphTimeStep(0.04*timeScale);
        logger(INFO,"Oomph time step %", getOomphTimeStep());
        //set dissipation
        double criticalDissipation = 0.1 * sqrt(elasticModulus_ * density_ / lengthScale); //TW
        setDissipation(criticalDissipation);
        logger(INFO,"Adding dissipation %",criticalDissipation);
        //setBodyForceAsGravity();
        setNewtonSolverTolerance(3e-8);
        prepareForSolve();
        solveSteady();
        getCentralDeflection();
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
        species->setStiffnessAndRestitutionCoefficient(1000, 0.5, 2.0*mass);
        logger(INFO,"Species: %", *species);
        // add particle
        SphericalParticle p(species);
        p.setRadius(radius);
        p.setPosition({getXCenter(),getYCenter(),getZMax()+p.getRadius()});
        //p.setPosition({getXCenter(),getYMax()+p.getRadius(),getZCenter()});
        p.setVelocity({0,0,-1});
        auto particle = particleHandler.copyAndAddObject(p);
        logger(INFO,"Particle: %", *particle);
        // set time step
        double tc = species->getCollisionTime(species->getMassFromRadius(p.getRadius()));
        setTimeStep(0.02*tc);
        logger(INFO,"Mercury time step %", getTimeStep());
        // set output file properties
        setParticlesWriteVTK(true);
        setWallsWriteVTK(true);
        setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        restartFile.writeFirstAndLastTimeStep();
        // add gravity
        setGravity({0,0,-0});

    }

    void actionsBeforeOomphTimeStep() override {
        static std::ofstream out(getName()+".out");
        out << getTime() << ' ' << getCentralDeflection() << ' ' << particleHandler.getLastObject()->getPosition().Z << '\n';
    }

    /**
     * Outputs deflection at midpoint
     */
    double getCentralDeflection() {
        std::array<double,3> min;
        std::array<double,3> max;
        getDomainSize(min, max);
        
        Vector<double> xi(3);
        xi[0] = 0.5*(max[0]+min[0]);
        xi[1] = 0.5*(max[1]+min[1]);
        xi[2] = 0.5*(max[2]+min[2]);
        double deflection = getDeflection(xi,2);
        logger(INFO, "Beam deflection at center (% % %) is %",
               xi[0], xi[1],xi[2], deflection);

        double mass, elasticEnergy, kineticEnergy;
        Vector<double> linearMomentum(3), angularMomentum(3);
        getMassMomentumEnergy(mass, linearMomentum, angularMomentum, elasticEnergy, kineticEnergy);
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
