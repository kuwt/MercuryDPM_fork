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
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <navier_stokes/refineable_navier_stokes_elements.h>
#include "Oomph/VCoupling/VCoupledFluidProblem.h"

/**
 * Define a particle-fluid coupled problem
 */
class CoupledSettlingParticle : public VCoupledFluidProblem<RefineableQCrouzeixRaviartElement<3>> //FIXME This does not contain AJ equations but Navier-Stokes
{
public:
    void setupOomph()
    {
        setRefineableFluidCubicMesh(10, 10, 10, 0.0, 1.0, 0.0, 1.0, 0.0, 5.0);
        pinBC(0);
        pinBC(1);
        pinBC(2);
        pinBC(3);
        pinBC(4);
        setBC(0, 2, 2.40);
        
        setReynoldsNumber(1.0);
        setReynoldsStrouhalNumber(1.0);
        setReynolds_InverseFroudeNumber(1.0);
        
        prepareForSolve();
    }
    
    void setupMercury()
    {
        //set domain
        setDomain({0.0, 0.0, 0.0}, {1.0, 1.0, 5.0});
    
        // add elastic species
        const double radius = 0.02;
    
        // add elastic species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(2.0e4);
        species->setStiffness(1.0e8*2*radius);
        SphericalParticle p(species);
        p.setRadius(radius);
        const double mass = species->getMassFromRadius(radius);
        double overlap = mass*9.81/species->getStiffness();
        p.setPosition({0.55, 0.55, 3.751});
        p.setVelocity({0.0,0.0,0.0});
        auto particle = particleHandler.copyAndAddObject(p);
        
        p.setRadius(1.5*radius);
        p.setPosition({0.45, 0.45, 3.751});
        p.setVelocity({0.0,0.0,0.0});
        particleHandler.copyAndAddObject(p);

        // set time step
        double tc = species->getCollisionTime(species->getMassFromRadius(p.getRadius()));
        setTimeStep(0.05*tc);
        
        // set oomph and mercury time step
        setOomphTimeStep(1000.0*getTimeStep());
        logger(INFO,"Mercury time step %", getTimeStep());
        logger(INFO,"Oomph time step %", getOomphTimeStep());
        
        // set output file properties
        setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::ONE_FILE);
        setSaveCount(100);
        restartFile.setFileType(FileType::ONE_FILE);
        restartFile.writeFirstAndLastTimeStep();
        setParticlesWriteVTK(true);
        // add gravity
        setGravity({0,0,-9.81});
    }
    
    void actionsBeforeOomphTimeStep() override
    {
        updateFluidVolumeFractions();
    };
    
    void actionsBeforeTimeStep() override
    {
        updateFluidVolumeFractions();
    }

protected:

    
};

int main()
{
    // Creating problem class and setting some properties
    CoupledSettlingParticle problem;
    problem.setName("CoupledSettlingParticle");
    
    problem.setTimeMax(3.0*problem.getTimeStep());
    logger(INFO,"timeMax: %, nTimeSteps %", problem.getTimeMax(), problem.getTimeMax()/problem.getTimeStep());
    problem.linear_solver_pt()->disable_doc_time();
    
    // Set fluid properties
    problem.setFluidDensity(1000.0);
    problem.setFluidDynamicViscosity(8.91e-4);
    problem.setFluidKinematicViscosity(8.91e-4);
    
    // Initialize Oomph and Mercury problems
    problem.setupOomph();
    problem.setupMercury();
    
    // Interactions (Not 4-way coupled due to NS-fluid-element instead of AJ-fluid-element)
    problem.setFractionDragForce(1.0);
    problem.updateFluidVolumeFractions();

    // Setting output options
    DocInfo docInfo;
    problem.doc_solution(docInfo);
    docInfo.number()++;
    
    // Solve the  problem
    problem.steady_newton_solve(); // Perform 1 steady newton solve before simulation starts
    problem.doc_solution(docInfo);
    docInfo.number()++;
    // Time-loop
    problem.solveParticleFluidCoupling();
    
    // Output after last timestep
    problem.doc_solution(docInfo);
    problem.doc_paraview(docInfo);
    
    return 0;
}