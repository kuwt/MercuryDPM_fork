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

#include "SolidProblem.h"

class SolidBeam : public SolidProblem
{
public:
    
    /**
     * https://autofem.com/examples/deflection_of_a_plate_under_theg.html
     */
    void getBeamDeflection()
    {
        
        std::array<double,3> min;
        std::array<double,3> max;
        getDomainSize(min, max);
        
        Vector<double> xi(3);
        xi[0] = max[0];
        xi[1] = 0.5*(max[1]+min[1]);
        xi[2] = 0.5*(max[2]+min[2]);
        double deflection = getDeflection(xi,2);
        logger(INFO, "Beam deflection at right end (% % %) is %",
               xi[0], xi[1],xi[2], deflection);
        
        double length = max[0]-min[0];
        double height = max[2]-min[2];
        logger(INFO, "Deflection should be -3*rho*g*L^4/(2*E*h^2) = % (in practice: %)",
               -3. * 9.8 * density_ * pow(length, 4) / ( 2. * elasticModulus_ * pow(height, 2) ), -0.00131642);
    }
};

/**
 * Measure bag deformation loaded by body force
 */
int main()
{
    // Solve the problem
    SolidBeam problem;
    problem.setName("SolidBeam");
    problem.setElasticModulus(1e8);
    problem.setDensity(2500);
    problem.setSolidCubicMesh(20, 2, 2, 0, 0.2, 0, 0.02, 0, 0.02);
    problem.setIsPinned([](SolidNode* n, unsigned d) {
        return n->is_on_boundary(SolidProblem::Boundary::X_MIN);
    });
    problem.setGravityAsBodyForce();
    problem.set_newton_solver_tolerance(3e-8);
    problem.prepareForSolve();
    problem.writeToVTK();
    problem.solveSteady();
    problem.getBeamDeflection();
    return 0;
}
