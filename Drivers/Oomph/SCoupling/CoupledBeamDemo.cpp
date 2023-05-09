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

//! [CoupledBeamDemo]
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Oomph/SCoupling/SCoupledSolidProblem.h"

int main()
{
    SCoupledSolidProblem<RefineableQDPVDElement<3, 2>> problem;
    //set name
    problem.setName("CoupledBeamUnitTest");

    //setup mercury
    problem.setDomain({0,0,0},{0.16,0.04,0.01});
    // add species
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    species->setDensity(1000);
    species->setStiffness(5.2);
    // add particle
    SphericalParticle p(species);
    p.setRadius(5e-3);
    p.setPosition({problem.getXCenter(),problem.getYCenter(),problem.getZMax()+p.getRadius()});
    p.setVelocity({0,0,-0.05});
    auto particle = problem.particleHandler.copyAndAddObject(p);
    // set time step
    double tc = species->getCollisionTime(species->getMassFromRadius(p.getRadius()));
    problem.setTimeStep(0.02*tc);

    // setup oomph
    problem.setElasticModulus(1e6);
    problem.setDensity(2500);
    problem.setSolidCubicMesh(7, 5, 5, 0, 0.16, 0, 0.04, 0, 0.01);
    problem.pinBoundaries({problem.Boundary::X_MIN, problem.Boundary::X_MAX});
    problem.setNewtonSolverTolerance(3e-8);
    problem.prepareForSolve();
    problem.coupleBoundary(problem.Boundary::Z_MAX);
    // sync time steps
    problem.setOomphTimeStep(problem.getTimeStep());

    // setup time
    problem.setTimeMax(100*problem.getTimeStep());

    // Solve the  problem and save mesh
    problem.solveSurfaceCoupling();
    problem.saveSolidMesh();
    return 0;
}
//! [CoupledBeamDemo]
