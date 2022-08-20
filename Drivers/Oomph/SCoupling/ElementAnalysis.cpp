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
#include "Coupling/SurfaceCoupledSolidProblem.h"

/**
 * Define a coupled problem
 */
class ElementAnalysis : public SurfaceCoupledSolidProblem<RefineableQDPVDElement<3, 2>>
{
public:
    ElementAnalysis () {
        //set name
        setName("ElementAnalysis");
        #ifdef OOMPH_HAS_MUMPS
        linear_solver_pt()=new MumpsSolver;
        #endif
        //remove existing output files
        removeOldFiles();

        // setup oomph
        setElasticModulus(1e6);
        setDensity(2500);
        setSolidCubicMesh(3, 3, 3, -1, 1, -1, 1, -1, 1);
        pinBoundaries({Boundary::X_MIN, Boundary::X_MAX});
        prepareForSolve();
        coupleBoundary(Boundary::Z_MAX);

        // number of elements
        logger(INFO, "\nElements: %", solid_mesh_pt()->nelement());
        std::cout << "Center Positions:\n";
        for (int e = 0; e < solid_mesh_pt()->nelement(); ++e) {
            ELEMENT* e_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(e));
            std::cout << " " << e << ": " << getXiCenter(e_pt) << std::endl;
        }

        ELEMENT* e_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(22));
//        // Set up memory for the shape/test functions
//        Shape psi(e_pt->nnode());
//        Vector<double> s(3, 0.0);
//        // Get shape/test fcts
//        e_pt->shape(s, psi);
//        GeomObject* geom_obj_pt = 0;
//        e_pt->locate_zeta(x, geom_obj_pt, s);

        logger(INFO, "Nodes per element: %", e_pt->nnode());
        std::cout << "Nodal positions element 22:\n";
        for (int n = 0; n < e_pt->nnode(); ++n) {
            SolidNode* n_pt = dynamic_cast<SolidNode*>(e_pt->node_pt(n));
            std::cout << " " << n << ": " << n_pt->x(0) << ' ' << n_pt->x(1) << ' ' << n_pt->x(2) << std::endl;
        }

        // Loop over the integration points
        logger(INFO, "Integration points per element: %", e_pt->integral_pt()->nweight());
        std::cout << "Local and global coordinates, weight, Jacobian and shape functions at integration points:\n";
        unsigned dim = e_pt->nodal_dimension();
        for (unsigned ipt = 0; ipt < e_pt->integral_pt()->nweight(); ipt++) {
            Vector<double> s(dim);
            // Assign the values of s
            for (unsigned i = 0; i < dim; ++i) {
                s[i] = e_pt->integral_pt()->knot(ipt, i);
            }
            //set shape and its derivatives
            Shape psi(e_pt->nnode(), e_pt->nnodal_position_type());
            DShape dpsidxi(e_pt->nnode(), e_pt->nnodal_position_type(), dim);
            double J = e_pt->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);
            // weight
            double w = e_pt->integral_pt()->weight(ipt);
            //interpolate
            Vector<double> interpolated_x(dim, 0.0);
            Vector<double> interpolated_xi(dim, 0.0);
            std::stringstream psi_ss;
            for (int l = 0; l < e_pt->nnode(); ++l) { // test functions
                for (int k = 0; k < e_pt->nnodal_position_type(); ++k) { // eqn_number
                    for (int i = 0; i < dim; ++i) {
                        //local_eqn = position_local_eqn_at_node(k, i);
                        interpolated_xi[i] += e_pt->lagrangian_position_gen(l, k, i) * psi(l,k);
                        interpolated_x[i] += e_pt->nodal_position_gen(l, k, i) * psi(l, k);
                    }
                    psi_ss << psi(l, k) << ' ';
                }
            }
            std::cout << "ipt " << ipt << " w " << w << " J " << J << " s " << s << " xi " << interpolated_xi << " psi " << psi_ss.str() << std::endl;
        }
    }

    static Vector<double> getXiCenter(ELEMENT* e_pt) {
        Vector<double> center(3);
        for (int n = 0; n < e_pt->nnode(); ++n) {
            SolidNode* n_pt = dynamic_cast<SolidNode*>(e_pt->node_pt(n));
            for (unsigned d = 0; d < n_pt->ndim(); ++d) {
                center[d] += n_pt->xi(d) / e_pt->nnode();
            }
        }
        return center;
    }
};

int main()
{
    ElementAnalysis problem;
    return 0;
}
