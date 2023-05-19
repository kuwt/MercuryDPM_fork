
//////////////////////////////////////////////////////////////
//  abandoned preliminary version; not maintained
//////////////////////////////////////////////////////////////


//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
// Driver formesh gluing


// Generic oomph-lib headers
#include "generic.h"

// Solid mechanics
#include "solid.h"

// The mesh 
#include "meshes/simple_cubic_mesh.h"

using namespace std;

using namespace oomph;


// Temporary inclusion; to be tiedied up into proper header
#include "glued_mesh_stuff.h"


//================================================================
/// Global variables
//================================================================
namespace Global_Physical_Variables {
    /// Pointer to strain energy function
    StrainEnergyFunction *Strain_energy_function_pt;

    /// Pointer to constitutive law
    ConstitutiveLaw *Constitutive_law_pt;

    /// Elastic modulus
    double E = 1.0;

    /// Poisson's ratio
    double Nu = 0.3;

    /// "Mooney Rivlin" coefficient for generalised Mooney Rivlin law
    double C1 = 1.3;

    /// Uniform pressure
    double P = 0.0;

    /// Constant pressure load
    void constant_pressure(const Vector<double> &xi, const Vector<double> &x,
                           const Vector<double> &n, Vector<double> &traction) {
        unsigned dim = traction.size();
        for (unsigned i = 0; i < dim; i++) {
            traction[i] = -P * n[i];
        }

    } // end of pressure load

}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



//====================================================================== 
/// Deformation of elastic pouch
//====================================================================== 
template<class ELEMENT>
class SheetGlueProblem : public Problem {

public:

    /// Constructor:
    SheetGlueProblem();

    /// Run simulation.
    void run(const std::string &dirname);

    /// Doc the solution
    void doc_solution(DocInfo &doc_info);

    /// Update function (empty)
    void actions_after_newton_solve() {}


    /// Create traction elements
    void create_traction_elements(SolidMesh *traction_mesh_pt,
                                  SolidMesh *solid_mesh_pt,
                                  const unsigned &b) {
        unsigned n_element = solid_mesh_pt->nboundary_element(b);
        for (unsigned e = 0; e < n_element; e++) {
            // The element itself:
            FiniteElement *fe_pt = solid_mesh_pt->boundary_element_pt(b, e);

            // Find the index of the face of element e along boundary b
            int face_index = solid_mesh_pt->face_index_at_boundary(b, e);


            // Create new element
            traction_mesh_pt->add_element_pt(
                    new SolidTractionElement<ELEMENT>
                            (fe_pt, face_index));

        }

        // Complete build process for SolidTractionElements
        n_element = traction_mesh_pt->nelement();
        for (unsigned i = 0; i < n_element; i++) {

            //Cast to a solid traction element
            SolidTractionElement<ELEMENT> *el_pt =
                    dynamic_cast<SolidTractionElement<ELEMENT> *>
                    (traction_mesh_pt->element_pt(i));

            //Set the traction function
            el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
        }

    }

    /// Update before solve: Empty
    void actions_before_newton_solve() {}

private:

    /// Pointers to solid meshes
    Vector<SolidMesh *> Solid_mesh_pt;

    /// Pointers to meshes of traction elements
    Vector<SolidMesh *> Traction_mesh_pt;

    /// Glued mesh
    GluedSolidMesh *Glued_mesh_pt;

};

//====================================================================== 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
SheetGlueProblem<ELEMENT>::SheetGlueProblem() {

    bool glue_it = true;

    Problem::Max_newton_iterations = 20;

    // make space for meshes
    unsigned n_mesh = 2;
    Solid_mesh_pt.resize(n_mesh);
    Traction_mesh_pt.resize(n_mesh);

    // How thick is the sheet?
    double sheet_thickness = 0.1;

    // Sheet dimensions
    double x_min = -sheet_thickness / 2.0;
    double x_max = sheet_thickness / 2.0;
    double y_min = -1.0;
    double y_max = 1.0;
    double z_min = -1.0;
    double z_max = 1.0;

    // Elements
    unsigned n_x = 2;
    unsigned n_y = 5;
    unsigned n_z = 5;

    // Width of glued layer:
    double glue_layer_width_y = 2.0 / double(n_y);
    double glue_layer_width_z = 2.0 / double(n_z);


    // Create bulk meshes
    for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++) {
        //Now create the mesh
        Solid_mesh_pt[i_mesh] = new SolidCubicMesh<ELEMENT>(n_x, n_y, n_z,
                                                            x_min, x_max,
                                                            y_min, y_max,
                                                            z_min, z_max);

        //Output boundaries
        std::string filename = "mesh_boundaries" + to_string(i_mesh) + ".dat";
        Solid_mesh_pt[i_mesh]->output_boundaries(filename);

        // Solid mesh is first sub-mesh
        if (!glue_it) {
            add_sub_mesh(Solid_mesh_pt[i_mesh]);
        }

        //Loop over the elements in the mesh to set parameters/function pointers
        unsigned n_element = Solid_mesh_pt[i_mesh]->nelement();
        for (unsigned i = 0; i < n_element; i++) {
            //Cast to a solid element
            ELEMENT *el_pt =
                    dynamic_cast<ELEMENT *>(Solid_mesh_pt[i_mesh]->element_pt(i));

            // Set the constitutive law
            el_pt->constitutive_law_pt() =
                    Global_Physical_Variables::Constitutive_law_pt;
        }


        // Setup boundary conditions
        //--------------------------

        // hierher: use
        bool minimal_constraint = true;

        // Suppress rigid body modes by pinning six displacements
        //--------------------------------------------------------
        // at lower-left, etc. nodes on
        //------------------------------
        // const y face at the glued end of the first sheet
        //-------------------------------------------------
        if (i_mesh == 0) {

            unsigned b = 4;
            SolidNode *ll_node_pt = Solid_mesh_pt[i_mesh]->boundary_node_pt(b, 0);
            SolidNode *ur_node_pt = Solid_mesh_pt[i_mesh]->boundary_node_pt(b, 0);
            SolidNode *ul_node_pt = Solid_mesh_pt[i_mesh]->boundary_node_pt(b, 0);
            SolidNode *lr_node_pt = Solid_mesh_pt[i_mesh]->boundary_node_pt(b, 0);

            unsigned n_node = Solid_mesh_pt[i_mesh]->nboundary_node(b);
            for (unsigned n = 0; n < n_node; n++) {
                SolidNode *nod_pt = Solid_mesh_pt[i_mesh]->boundary_node_pt(b, n);
                double y = nod_pt->x(1);
                double z = nod_pt->x(2);
                if ((y <= ll_node_pt->x(1)) && (z <= ll_node_pt->x(2))) {
                    ll_node_pt = nod_pt;
                }
                if ((y <= ul_node_pt->x(1)) && (z >= ul_node_pt->x(2))) {
                    ul_node_pt = nod_pt;
                }
                if ((y >= ur_node_pt->x(1)) && (z >= ur_node_pt->x(2))) {
                    ur_node_pt = nod_pt;
                }
                if ((y >= lr_node_pt->x(1)) && (z <= lr_node_pt->x(2))) {
                    lr_node_pt = nod_pt;
                }
            }

            // Note sharp deformation near ll node because it has to balance
            // the entire z force exerted downwards by the pressure acting on the
            // increasingly inflated pouch

            // Suppress rigid body motion
            ll_node_pt->pin_position(0);
            ll_node_pt->pin_position(1);
            ll_node_pt->pin_position(2);

            // Suppress rotation about x and y axes
            ul_node_pt->pin_position(0);
            ul_node_pt->pin_position(1);

            // suppress rotation about z axis
            lr_node_pt->pin_position(0);

            oomph_info << "Pinning at: \n"
                       << ll_node_pt->x(0) << " "
                       << ll_node_pt->x(1) << " "
                       << ll_node_pt->x(2) << "\n"
                       << ul_node_pt->x(0) << " "
                       << ul_node_pt->x(1) << " "
                       << ul_node_pt->x(2) << "\n"
                       << lr_node_pt->x(0) << " "
                       << lr_node_pt->x(1) << " "
                       << lr_node_pt->x(2) << "\n";

        }

        //######################################################################


        // // y=0
        // //----
        // {
        //  unsigned b=1;
        //  unsigned n_node = Solid_mesh_pt[i_mesh]->nboundary_node(b);
        //  for(unsigned n=0;n<n_node;n++)
        //   {
        //    //Pin all nodes
        //    for(unsigned i=0;i<3;i++)
        //     {
        //      Solid_mesh_pt[i_mesh]->boundary_node_pt(b,n)->pin_position(i);
        //     }
        //   }
        // }
        // // y=1
        // //----
        // {
        //  unsigned b=3;
        //  unsigned n_node = Solid_mesh_pt[i_mesh]->nboundary_node(b);
        //  for(unsigned n=0;n<n_node;n++)
        //   {
        //    //Pin all nodes
        //    for(unsigned i=0;i<3;i++)
        //     {
        //      Solid_mesh_pt[i_mesh]->boundary_node_pt(b,n)->pin_position(i);
        //     }
        //   }
        // }

        //######################################################################

        // Pin the redundant solid pressures
        PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
                Solid_mesh_pt[i_mesh]->element_pt());
    }


    // Shift second sheet:
    unsigned nnode = Solid_mesh_pt[1]->nnode();
    for (unsigned j = 0; j < nnode; j++) {
        Solid_mesh_pt[1]->node_pt(j)->x(0) += sheet_thickness;
    }
    Solid_mesh_pt[1]->set_lagrangian_nodal_coordinates();


    if (glue_it) {
        // Create glued mesh
        Glued_mesh_pt = new GluedSolidMesh(Solid_mesh_pt);


        // Select nodes to be glued:
        //--------------------------

        // We retain them on mesh 0 where the glued
        // nodes are on boundary 2
        {
            double tol = 1.0e-2;
            unsigned i_mesh = 0;
            unsigned b = 2;
            unsigned n_node = Solid_mesh_pt[i_mesh]->nboundary_node(b);
            oomph_info << "# of glue nodes: " << n_node << std::endl;
            Vector<Node *> glue_node_pt;
            glue_node_pt.reserve(n_node);
            for (unsigned n = 0; n < n_node; n++) {
                Node *potential_node_pt = Solid_mesh_pt[i_mesh]->boundary_node_pt(b, n);
                if (
                        (potential_node_pt->x(1) < (y_min + glue_layer_width_y + tol)) ||
                        (potential_node_pt->x(1) > (y_max - glue_layer_width_y - tol)) ||
                        (potential_node_pt->x(2) < (z_min + glue_layer_width_z + tol))) {
                    glue_node_pt.push_back(potential_node_pt);
                }
            }

            // Glue specified nodes to co-located nodes in consituent mesh 1
            unsigned i_mesh_replace = 1;
            Glued_mesh_pt->glue(glue_node_pt, i_mesh_replace);

            // Add to global mesh
            add_sub_mesh(Glued_mesh_pt);
        }
    } else {
        Glued_mesh_pt = 0;
    }


    // Create traction meshes (attach to elements in still existing
    //-------------------------------------------------------------
    // original meshes)
    //-----------------
    for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++) {
        // Make new mesh
        Traction_mesh_pt[i_mesh] = new SolidMesh;

        // Create elements on x=1
        unsigned b = 2;
        if (i_mesh == 1) {
            b = 4;
        }
        create_traction_elements(Traction_mesh_pt[i_mesh], Solid_mesh_pt[i_mesh], b);

        // Add it!
        add_sub_mesh(Traction_mesh_pt[i_mesh]);
    }


    // Build combined "global" mesh
    build_global_mesh();


    //Attach the boundary conditions to the mesh
    cout << assign_eqn_numbers() << std::endl;

}


//==================================================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void SheetGlueProblem<ELEMENT>::doc_solution(DocInfo &doc_info) {

    ofstream some_file;
    char filename[100];

    // Number of plot points
    unsigned npts = 5;

    // Output shape of deformed body
    sprintf(filename, "%s/sol1n%i.dat", doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    Solid_mesh_pt[0]->output(some_file, npts);
    Solid_mesh_pt[1]->output(some_file, npts);
    some_file.close();

    // Output traction
    //----------------
    sprintf(filename, "%s/traction%i.dat", doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    Traction_mesh_pt[0]->output(some_file, npts);
    Traction_mesh_pt[1]->output(some_file, npts);
    some_file.close();


    // Output glued mesh
    //------------------
    if (Glued_mesh_pt != 0) {
        sprintf(filename, "%s/glued_sol1n%i.dat", doc_info.directory().c_str(),
                doc_info.number());
        some_file.open(filename);
        Glued_mesh_pt->output(some_file, npts);
        some_file.close();
    }

}


//==================================================================
/// Run the problem
//==================================================================
template<class ELEMENT>
void SheetGlueProblem<ELEMENT>::run(const std::string &dirname) {

    // Output
    DocInfo doc_info;

    // Set output directory
    doc_info.set_directory(dirname);

    // Step number
    doc_info.number() = 0;

    // Doc initial configuration
    doc_solution(doc_info);
    doc_info.number()++;


    // Pressure increment
    double dp = 0.0005;
    Global_Physical_Variables::P = 0.0; // dp;

    //Parameter incrementation
    unsigned nstep = 10;
    if (CommandLineArgs::Argc != 1) {
        std::cout << "Validation -- only doing one step" << std::endl;
        nstep = 1;
    }

    for (unsigned i = 0; i < nstep; i++) {

        // Solve the problem with Newton's method
        newton_solve();

        // Doc solution
        doc_solution(doc_info);
        doc_info.number()++;

        //Increase the pressure
        Global_Physical_Variables::P += dp;
    }

}


//======================================================================
/// Driver for simple elastic problem
//======================================================================
int main(int argc, char **argv) {

    // Store command line arguments
    CommandLineArgs::setup(argc, argv);

    //Initialise physical parameters
    Global_Physical_Variables::E = 2.1;
    Global_Physical_Variables::Nu = 0.4;
    Global_Physical_Variables::C1 = 1.3;


    // Define a strain energy function: Generalised Mooney Rivlin
    Global_Physical_Variables::Strain_energy_function_pt =
            new GeneralisedMooneyRivlin(&Global_Physical_Variables::Nu,
                                        &Global_Physical_Variables::C1,
                                        &Global_Physical_Variables::E);

    // Define a constitutive law (based on strain energy function)
    Global_Physical_Variables::Constitutive_law_pt =
            new IsotropicStrainEnergyFunctionConstitutiveLaw(
                    Global_Physical_Variables::Strain_energy_function_pt);

    //Set up the problem with pure displacement formulation
    SheetGlueProblem<QPVDElement<3, 3> > problem;
    problem.run("./");


}





