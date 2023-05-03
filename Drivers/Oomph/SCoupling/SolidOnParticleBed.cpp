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
#include "Mercury3D.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "Oomph/SCoupling/SCoupledSolidProblem.h"

class CoupledProblem : public SCoupledSolidProblem<RefineableQDPVDElement<3, 2>>
{
public:
    
    CoupledProblem () {
        //set name
        setName("SolidOnParticleBed");
        //remove existing output files
        removeOldFiles();
        
        // setup steps
        setupOomph();
        setupMercury();
        setOomphTimeStep(getTimeStep());
        setTimeMax(2000*getTimeStep());
        setSaveCount(10);
        logger(INFO,"Time step %", getTimeStep());

        // Solve the  problem
        writeToVTK();
        solveSurfaceCoupling();
        saveSolidMesh();

        helpers::writeToFile(getName()+".gnu","set key autotitle columnheader\n"
                                              "p '" + getName() + ".def' u 1:2");
    }

    void setupOomph() {
        // set stiffness
        setElasticModulus(1e6);
        // set density
        setDensity(1000);
        // set solid dimensions
        unsigned nx=2, ny=2, nz=2;
        double xMax=5e-2, xMin=-xMax; // one litre of material, density 1000 -> 1 kg
        double yMax=5e-2, yMin=-yMax;
        double zMax=10e-2, zMin=0;
        setSolidCubicMesh(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax);
        // no pinned boundaries
        setIsPinned([](SolidNode* n, unsigned d) {
            return false;
        });
        // set dissipation
        setDissipation(0.0);
        // set gravity
        setBodyForceAsGravity();
        // set solver properties
        setNewtonSolverTolerance(3e-8);
        disable_info_in_newton_solve();
        linear_solver_pt()->disable_doc_time();
        // assemble
        prepareForSolve();
        // couple all boundaries
        coupleBoundary(Boundary::Z_MIN);

        double mass, elasticEnergy, kineticEnergy;
        Vector<double> com(3), linearMomentum(3), angularMomentum(3);
        getMassMomentumEnergy(mass, com, linearMomentum, angularMomentum, elasticEnergy, kineticEnergy);
        logger(INFO, "mass %, linearMomentum %, angularMomentum %, elasticEnergy %, totalEnergy %",
               mass, linearMomentum, angularMomentum, elasticEnergy, elasticEnergy+kineticEnergy);
    }

    void setupMercury() {
        //set domain for visualisation
        std::array<double, 3> min;
        std::array<double, 3> max;
        SolidProblem::getDomainSize(min, max);
        const double radius = 10e-3; //10mm radius
        min[2] = -2.5*radius;
        Mercury3D::setDomain(min, max);
        // add elastic species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(1000);
        species->setStiffness(9.8*250.0); //should result in 1mm overlap
        const double mass = species->getMassFromRadius(radius);
        double tc = species->getCollisionTime(mass);
        species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, 1, 1, mass);
        species->setSlidingFrictionCoefficient(0.2);
        // set time step
        setTimeStep(0.05*tc);
        // add particle
        SphericalParticle p(species);
        p.setRadius(radius);
        double x = 0.5*max[0], y = 0.8*max[1];
        p.setPosition({x, y, -radius});
        particleHandler.copyAndAddObject(p);
        p.setPosition({-x, y, -radius});
        particleHandler.copyAndAddObject(p);
        p.setPosition({x, -y, -radius});
        particleHandler.copyAndAddObject(p);
        p.setPosition({-x, -y, -radius});
        particleHandler.copyAndAddObject(p);
        // add bottom wall
        InfiniteWall w(species);
        w.set({0,0,-1},{0,0,-2.0*radius});
        wallHandler.copyAndAddObject(w);
        // set output file properties
        setParticlesWriteVTK(true);
        wallHandler.setWriteVTK(true);
        setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        restartFile.writeFirstAndLastTimeStep();
        dataFile.setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::ONE_FILE);
        // add gravity
        setGravity({0,0,-9.8});
    }

    void actionsAfterTimeStep () override {
        getSolidDeflection();
    }

    void getSolidDeflection() {
        std::array<double,3> min;
        std::array<double,3> max;
        getDomainSize(min, max);
        
        Vector<double> xi(3);
        xi[0] = 0.5*(max[0]+min[0]);
        xi[1] = 0.5*(max[1]+min[1]);
        xi[2] = min[2];
        double deflection = getDeflection(xi,2);

        static std::ofstream file(getName()+".def");
        file << getTime() << " " << deflection << std::endl;
    }
};

/**
 * Measure bag deformation loaded by body force
 */
int main()
{
    CoupledProblem problem;
    return 0;
}
