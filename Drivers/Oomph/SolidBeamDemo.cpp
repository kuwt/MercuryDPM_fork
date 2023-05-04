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

/**
 * Simulates a cubic rod (20mm x 2mm x 2mm) of density 2500 kg/m^3, elastic modulus 100 MPa, fixed on left side, bent by gravity
 */

//! [SolidBeamDemo]
#include "Oomph/SolidProblem.h"

int main()
{
    // Defines a SolidProblem of element type RefineableQDPVDElement<3,2>.
    SolidProblem<RefineableQDPVDElement<3, 2>> problem;
    // Sets up the problem
    problem.setName("SolidBeamUnitTest");
    problem.setElasticModulus(1e8);
    problem.setDensity(2500);
    problem.setSolidCubicMesh(20, 2, 2, 0, 0.2, 0, 0.02, 0, 0.02);
    problem.pinBoundary(problem.Boundary::X_MIN);
    problem.setBodyForceAsGravity();
    problem.setNewtonSolverTolerance(3e-8);
    problem.prepareForSolve();
    // Solve the problem
    problem.solveSteady();
    return 0;
}
//! [SolidBeamDemo]