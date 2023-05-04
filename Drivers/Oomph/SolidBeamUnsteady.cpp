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

#include "Oomph/SolidProblem.h"
#include "Math/Helpers.h"

/// Defines a SolidProblem of element type RefineableQDPVDElement<3,2>. Add functionality to write output
class Beam : public SolidProblem<RefineableQDPVDElement<3, 2>>
{
    /// output file stream
    std::ofstream out;

    /// Write header of output file
    void actionsBeforeSolve() override {
        helpers::writeToFile(getName()+".gnu", "set key autotitle columnheader\n"
           "p 'SolidBeamUnsteady.out' u 1:3 w l, '' u 1:4 w l, '' u 1:5 w l, '' u 1:($3+$4+$5) t 'totalEnergy' w l");
        out.open(getName()+".out");
        out << "time deflection elasticEnergy kineticEnergy gravEnergy\n";
    }

    /// Each time step, compute deflection, elastic, kinetic and gravitational energy, and write to output file
    void actionsBeforeOomphTimeStep() override {
        double mass, elasticEnergy, kineticEnergy;
        Vector<double> com(3), linearMomentum(3), angularMomentum(3);
        getMassMomentumEnergy(mass, com, linearMomentum, angularMomentum, elasticEnergy, kineticEnergy);
        static double comZ0 = com[2];
        double gravEnergy = 9.8*mass*(com[2]-comZ0);
        out << getOomphTime() << ' ' << getBeamDeflection() << ' ' << elasticEnergy << ' ' << kineticEnergy << ' ' << gravEnergy << std::endl;
        std::cout << getOomphTime() << ' ' << getBeamDeflection() << ' ' << elasticEnergy << ' ' << kineticEnergy << ' ' << gravEnergy << '\n';
    }

    /// Computes beam deflection
    double getBeamDeflection() const {
        std::array<double, 3> min, max;
        getDomainSize(min, max);

        Vector<double> xi(3);
        xi[0] = max[0];
        xi[1] = 0.5 * (max[1] + min[1]);
        xi[2] = 0.5 * (max[2] + min[2]);
        return getDeflection(xi, 2);
    }

};

/**
 * Simulates a cubic rod (20mm x 2mm x 2mm) of density 2500 kg/m^3, elastic modulus 100 MPa, fixed on left side, bent by gravity
 */
int main()
{
    // Solve the problem
    Beam problem;
    problem.setName("SolidBeamUnsteady");
    problem.setElasticModulus(1e8);
    problem.setDensity(2500);
    problem.setSolidCubicMesh(20, 2, 2, 0, 0.2, 0, 0.02, 0, 0.02);
    problem.pinBoundary(Beam::Boundary::X_MIN);
    problem.setBodyForceAsGravity();
    problem.setNewtonSolverTolerance(3e-8);
    problem.prepareForSolve();
    problem.linear_solver_pt()->disable_doc_time();
    //problem.disable_info_in_newton_solve();
    //set time scale of oscillation
    //https://vlab.amrita.edu/?sub=3&brch=175&sim=1080&cnt=1
    double b = 0.02, d = 0.02, l = 0.2;
    double inertia = b*d*d*d/12;
    double omega = 1.875*1.875*sqrt(problem.getElasticModulus()*inertia/(problem.getDensity()*b*d*std::pow(l,4)));
    double timeScale = 2*constants::pi/omega; //0.06
    // simulate two oscillations, check energy balance
    problem.solveUnsteady(2.0*timeScale, 0.04*timeScale, 10);
    return 0;
}