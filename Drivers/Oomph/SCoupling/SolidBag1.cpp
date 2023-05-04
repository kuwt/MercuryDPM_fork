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
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONrimIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONrimACT, SrimICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Oomph/SolidProblem.h"
#include "Math/Helpers.h"

class SolidBag : public SolidProblem<RefineableQDPVDElement<3, 2>>
{
    double gripperDisplacement = 0.0;
public:
    SolidBag()
    {
        //read in variables from parameter file
        double elasticModulus = 2e6;
        double poisson = 0.4;
        double density = 1000;
        double length = 0.16;
        double width = 0.13;
        double thickness = 1e-4;
        double resolution = 0.01;
        double timeMax = 0.75;
        double rim = 0.01;
        double gripperTraverse = 0.01;
        double gripperLength = 0.01;
        double gripperHeight = 0.01;
        double suctionForce = 1e3;

        setName("SolidBag1");
        removeOldFiles();
        setElasticModulus(elasticModulus);
        setPoissonRatio(poisson);
        setDensity(density);
        setSolidCubicMesh(width/resolution, 2, 2.0*length/resolution, -width / 2, width / 2, 0, thickness, -length, length);

        setIsPinned([length, gripperLength, gripperHeight, width, rim](SolidNode* n, unsigned d) {
            if (d == 0 || d == 2)
            { //pin gripper region in x- and z-direction
                return ( n->is_on_boundary(Boundary::X_MIN)
                    || n->is_on_boundary(Boundary::X_MAX) )
                    && fabs(n->xi(2) - ( length - gripperLength - 0.5*gripperHeight )) < 0.51*gripperHeight;
            }
            else /*if (d==1)*/ { //pin rim in y-direction
                return fabs(n->xi(0)) > width / 2 - 1.01 * rim || n->xi(2) < 1.01 * rim;
            }
        });
        
        // bend the plate into a bag
        bendBag();

        ////set time scale of oscillation
        ////https://vlab.amrita.edu/?sub=3&brch=175&sim=1080&cnt=1
        //double d = thickness, l = width/2.0;
        //double omega = sqrt(getElasticModulus()*std::pow(1.875*d,2)/(12*getDensity()*std::pow(l,4)));
        //double timeScale = 2*constants::pi/omega;
        setOomphTimeStep(0.01);
        logger(INFO, "Time step %", getOomphTimeStep());
        
        // add dissipation
        //addDissipation(getCriticalDissipation(width)); //todo: use dissipationScale?

        setSuctionAsBodyForce(suctionForce);

        setNewtonSolverTolerance(1e-6);
        linear_solver_pt()->disable_doc_time();
        //disable_info_in_newton_solve();
        prepareForSolve();

        solveSteady();

        // Solve the  problem
        gripperDisplacement = gripperTraverse / timeMax * getOomphTimeStep();
        solveUnsteady(timeMax, getOomphTimeStep(), 3);
    }

    /**
      * Squeeze the bag at the top from both sides by a given factor.
      * \param distance
      */
    void actionsBeforeOomphTimeStep() override
    {
        // number of nodes in mesh
        const unsigned n_node = solid_mesh_pt()->nnode();
        // for all nodes in mesh
        for (unsigned i = 0; i < n_node; i++)
        {
            SolidNode* n_pt = solid_mesh_pt()->node_pt(i);
            if (n_pt->position_is_pinned(0))
            {
                //factorMax = gripperDisplacement*dt/(0.5*W)
                //const double factor = factorMax * n_pt->xi(2) / (Global_Physical_Variables::length.Z);
                n_pt->x(0) += n_pt->x(0) > 0 ? -gripperDisplacement : gripperDisplacement;
            }
        }
    }
    
    void bendBag()
    {
        logger(INFO, "Bend bag");
        double gapThickness = 1e-6;
        // bend the plate into a bag
        Vector<double> xi(3);
        for (unsigned i = 0; i < solid_mesh_pt()->nnode(); i++)
        {
            SolidNode* n_pt = solid_mesh_pt()->node_pt(i);
            // shift all points in y by gapThickness/2
            n_pt->xi(1) += gapThickness / 2;
            n_pt->x(1) += gapThickness / 2;
            //rotate lower half (z<0) in yz-plane by 180 deg
            if (n_pt->xi(2) < 0)
            {
                n_pt->xi(2) *= -1;
                n_pt->xi(1) *= -1;
                n_pt->x(2) *= -1;
                n_pt->x(1) *= -1;
            }
            // center lowest row
            if (n_pt->xi(2) < gapThickness)
            {
                n_pt->xi(1) = 0;
                n_pt->x(1) = 0;
            }
        }
    }

    void setSuctionAsBodyForce(double bodyForce)
    {
        // define a static body force
        static double bodyForce_;
        bodyForce_ = bodyForce;
        logger(INFO, "Adding suction as body force, rho g = %", bodyForce);
        body_force_fct = [](const double& time, const Vector<double>& xi, Vector<double>& b) {
            b[0] = 0.0;
            b[1] = xi[1] >= 0 ? bodyForce_ : -bodyForce_;
            b[2] = 0.0;
        };
    }
};

int main() {
    SolidBag problem;
    return 0;
}
