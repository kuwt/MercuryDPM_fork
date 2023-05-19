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
namespace Global_Physical_Variables
{
    /// Pointer to strain energy function
    StrainEnergyFunction* Strain_energy_function_pt;

    /// Pointer to constitutive law
    ConstitutiveLaw* Constitutive_law_pt;

    /// Elastic modulus
    double E=1.0;

    /// Poisson's ratio
    double Nu=0.3;

    /// "Mooney Rivlin" coefficient for generalised Mooney Rivlin law
    double C1=1.3;

    /// Uniform pressure
    double P = 0.0;

    /// Constant pressure load
    void constant_pressure(const Vector<double> &xi,const Vector<double> &x,
                           const Vector<double> &n, Vector<double> &traction)
    {
        unsigned dim = traction.size();
        for(unsigned i=0;i<dim;i++)
        {
            traction[i] = -P*n[i];
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
class SheetGlueProblem : public Problem
{

public:

    /// Constructor:
    SheetGlueProblem();

    /// Run simulation.
    void run(const std::string &dirname);

    /// Doc the solution
    void doc_solution(DocInfo& doc_info);

    /// Update function (empty)
    void actions_after_newton_solve() {}

    /// Create traction elements
    void create_traction_elements(SolidMesh* traction_mesh_pt,
                                  SolidMesh* solid_mesh_pt,
                                  const unsigned& b)
    {
        unsigned n_element = solid_mesh_pt->nboundary_element(b);
        for (unsigned e=0;e<n_element;e++)
        {
            // The element itself:
            FiniteElement* fe_pt = solid_mesh_pt->boundary_element_pt(b,e);

            // Find the index of the face of element e along boundary b
            int face_index = solid_mesh_pt->face_index_at_boundary(b,e);

            // Create new element
            traction_mesh_pt->add_element_pt(
                    new SolidTractionElement<ELEMENT>
                            (fe_pt,face_index));

        }

        // Complete build process for SolidTractionElements
        n_element=traction_mesh_pt->nelement();
        for(unsigned i=0;i<n_element;i++)
        {

            //Cast to a solid traction element
            SolidTractionElement<ELEMENT> *el_pt =
                    dynamic_cast<SolidTractionElement<ELEMENT>*>
                    (traction_mesh_pt->element_pt(i));

            //Set the traction function
            el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
        }

    }


    /// Update before solve: Empty
    void actions_before_newton_solve(){}

private:



    /// Helper function to get unique nodes of central sheets
    /// where they meet the sidewalls. This happens after the colocated
    /// nodes have been glued together, so visiting them through the
    /// boundary lookup schemes of the glued mesh's constituent meshes
    /// produces duplicates.
    void get_unique_nodes_that_meet_sidewalls(const bool& on_left,
                                              std::set<SolidNode*>& unique_nodes)
    {

        // Wipe
        unique_nodes.clear();

        // On left (min. y) or right (max. y)?
        unsigned b=1;
        if (!on_left)
        {
            b=3;
        }

        // Populate set
        for (unsigned i_mesh=0;i_mesh<4;i_mesh++)
        {
            unsigned n_node = Solid_mesh_pt[i_mesh]->nboundary_node(b);
            for(unsigned n=0;n<n_node;n++)
            {
                unique_nodes.insert(Solid_mesh_pt[i_mesh]->boundary_node_pt(b,n));
            }
        }

        // oomph_info << "Unique nodes on ";
        // if (on_left)
        //  {
        //   oomph_info << "left : ";
        //  }
        // else
        //  {
        //   oomph_info << "right: ";
        //  }
        // oomph_info << unique_nodes.size() << std::endl;

    }

    /// Pointers to solid meshes
    Vector<SolidMesh*> Solid_mesh_pt;

    /// Pointers to meshes of traction elements
    Vector<SolidMesh*> Traction_mesh_pt;

    /// Glued mesh
    GluedSolidMesh* Glued_mesh_pt;

};

//====================================================================== 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
SheetGlueProblem<ELEMENT>::SheetGlueProblem()
{


    // Make space for meshes that make up the actual sheets
    unsigned n_mesh=4;
    Solid_mesh_pt.resize(n_mesh);
    Traction_mesh_pt.resize(4);

    // How thick is the sheet?
    double sheet_thickness=0.1;

    // Overall pouch dimensions (correct for first sheet)
    double x_min_first_sheet=0.0;
    double x_max_first_sheet=sheet_thickness;
    double y_min_first_sheet=0.0;
    double y_max_first_sheet=1.0;
    double z_min_first_sheet=0.0;
    double z_max_first_sheet=1.0;

    // nz element in bits at the bottom (at least three)
    unsigned n_z_in_bottom=3;

    // Thickness: at least two
    unsigned n_x_first_sheet = 2;

    // Width: Whatever; at lest three
    unsigned n_y_first_sheet = 5;

    // Overall height (mor than bottom bits)
    unsigned n_z_first_sheet = 5;

    // height of folds at bottom
    double fold_height=double(n_z_in_bottom)/
                       double(n_z_first_sheet)*z_max_first_sheet;

    // Create bulk meshes
    for (unsigned i_mesh=0;i_mesh<n_mesh;i_mesh++)
    {
        double x_min=x_min_first_sheet;
        double x_max=x_max_first_sheet;
        double y_min=y_min_first_sheet;
        double y_max=y_max_first_sheet;
        double z_min=z_min_first_sheet;
        double z_max=z_max_first_sheet;

        unsigned n_x=n_x_first_sheet;
        unsigned n_y=n_y_first_sheet;
        unsigned n_z=n_z_first_sheet;

        switch (i_mesh)
        {
            // Outer bits
            case 0:
            case 3:
                break;

                // Shorter inner bits making up the fold
            case 1:
            case 2:

                // Update
                x_min=x_min_first_sheet;
                x_max=x_max_first_sheet;
                y_min=y_min_first_sheet;
                y_max=y_max_first_sheet;
                z_min=z_min_first_sheet;
                z_max=z_min_first_sheet+fold_height;
                n_z=n_z_in_bottom;

                break;


            default:
                oomph_info << "Never get here" << std::endl;
                abort(); // hierher proper error, please
        }


        //Now create the mesh
        Solid_mesh_pt[i_mesh] = new SolidCubicMesh<ELEMENT>(n_x,n_y,n_z,
                                                            x_min,x_max,
                                                            y_min,y_max,
                                                            z_min,z_max);



        // Pin the redundant solid pressures
        PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
                Solid_mesh_pt[i_mesh]->element_pt());



        // Shift other sheets:
        if (i_mesh>0)
        {
            unsigned nnode=Solid_mesh_pt[i_mesh]->nnode();
            for (unsigned j=0;j<nnode;j++)
            {
                Solid_mesh_pt[i_mesh]->node_pt(j)->x(0)+=double(i_mesh)*sheet_thickness;
            }
            Solid_mesh_pt[i_mesh]->set_lagrangian_nodal_coordinates();
        }


        // Output boundaries -- use these to identify boundaries where
        // BCs (pinning or traction) have to be applied
        std::string filename="mesh_boundaries"+to_string(i_mesh)+".dat";
        Solid_mesh_pt[i_mesh]->output_boundaries(filename);

        filename="mesh"+to_string(i_mesh)+".dat";
        Solid_mesh_pt[i_mesh]->output(filename);

    }


    // Create sidewall meshes
    //-----------------------
    {
        double x_min=x_min_first_sheet;
        double x_max=4*sheet_thickness;
        double y_min=-sheet_thickness;
        double y_max=y_min_first_sheet;
        double z_min=z_min_first_sheet;
        double z_max=z_max_first_sheet;

        unsigned n_x=4*n_x_first_sheet;
        unsigned n_y=n_x_first_sheet;
        unsigned n_z=n_z_first_sheet;

        //Now create the mesh and add to collection
        Solid_mesh_pt.push_back(new SolidCubicMesh<ELEMENT>(n_x,n_y,n_z,
                                                            x_min,x_max,
                                                            y_min,y_max,
                                                            z_min,z_max));

        // hierher Manual counting. Bad!
        Solid_mesh_pt[4]->output("sidewall_mesh1.dat");


        //Now create the mesh and to collection
        Solid_mesh_pt.push_back(new SolidCubicMesh<ELEMENT>(n_x,n_y,n_z,
                                                            x_min,x_max,
                                                            y_min,y_max,
                                                            z_min,z_max));

        // Shift it hierher Manual counting. Bad!
        unsigned nnode=Solid_mesh_pt[5]->nnode();
        for (unsigned j=0;j<nnode;j++)
        {
            Solid_mesh_pt[5]->node_pt(j)->x(1)+=(y_max_first_sheet-y_min_first_sheet)+
                                                sheet_thickness;
        }
        Solid_mesh_pt[5]->set_lagrangian_nodal_coordinates();

        Solid_mesh_pt[5]->output("sidewall_mesh2.dat");

    }


    // finish off
    //-----------
    unsigned nmesh=Solid_mesh_pt.size();
    for (unsigned i_mesh=0;i_mesh<nmesh;i_mesh++)
    {
        //Loop over the elements in the mesh to set parameters/function pointers
        unsigned  n_element=Solid_mesh_pt[i_mesh]->nelement();
        for(unsigned i=0;i<n_element;i++)
        {
            //Cast to a solid element
            ELEMENT *el_pt =
                    dynamic_cast<ELEMENT*>(Solid_mesh_pt[i_mesh]->element_pt(i));

            // Set the constitutive law
            el_pt->constitutive_law_pt() =
                    Global_Physical_Variables::Constitutive_law_pt;
        }
    }

    // Create glued mesh
    Glued_mesh_pt = new GluedSolidMesh(Solid_mesh_pt);
    Glued_mesh_pt->output("glued_mesh.dat");

    // Select nodes to be glued:
    //--------------------------
    double tol=1.0e-2;

    // First sheet glued to second sheet
    //----------------------------------
    {

        // Width of glued layer:
        double glue_layer_width_z=
                (z_max_first_sheet-z_min_first_sheet)/double(n_z_first_sheet);

        // We retain them on mesh 1 where the glued
        // nodes are on boundary 4
        unsigned i_mesh=1;
        unsigned b=4;
        unsigned n_node = Solid_mesh_pt[i_mesh]->nboundary_node(b);
        Vector<Node*> glue_node_pt;
        glue_node_pt.reserve(n_node);
        std::set<Node*> glue_node_set_pt;
        for(unsigned n=0;n<n_node;n++)
        {
            Node* potential_node_pt=Solid_mesh_pt[i_mesh]->boundary_node_pt(b,n);

            bool in_glued_strip=(
                    (potential_node_pt->x(2)<(z_min_first_sheet+glue_layer_width_z+tol))||
                    (potential_node_pt->x(2)>(z_max_first_sheet-glue_layer_width_z-tol))  );

            bool on_edge_to_be_glued_to_sidewall=(
                    (potential_node_pt->x(1)<(y_min_first_sheet+tol))||
                    (potential_node_pt->x(1)>(y_max_first_sheet-tol))  );

            // Use set to avoid duplicates
            if (in_glued_strip||on_edge_to_be_glued_to_sidewall)
            {
                glue_node_set_pt.insert(potential_node_pt);
            }
        }

        // Transfer to vector
        for (std::set<Node*>::iterator it=glue_node_set_pt.begin();
             it!=glue_node_set_pt.end();it++)
        {
            glue_node_pt.push_back((*it));
        }

        // Glue specified nodes to co-located nodes in consituent mesh 0
        unsigned i_mesh_replace=0;
        Glued_mesh_pt->glue(glue_node_pt,i_mesh_replace);
    }

    // Fourth sheet glued to third sheet
    //----------------------------------
    {

        // Width of glued layer:
        double glue_layer_width_z=
                (z_max_first_sheet-z_min_first_sheet)/double(n_z_first_sheet);

        // We retain them on mesh 2 where the glued
        // nodes are on boundary 2
        unsigned i_mesh=2;
        unsigned b=2;
        unsigned n_node = Solid_mesh_pt[i_mesh]->nboundary_node(b);
        Vector<Node*> glue_node_pt;
        std::set<Node*>  glue_node_set_pt;
        glue_node_pt.reserve(n_node);
        for(unsigned n=0;n<n_node;n++)
        {
            Node* potential_node_pt=Solid_mesh_pt[i_mesh]->boundary_node_pt(b,n);

            bool in_glued_strip=(
                    (potential_node_pt->x(2)<(z_min_first_sheet+glue_layer_width_z+tol))||
                    (potential_node_pt->x(2)>(z_max_first_sheet-glue_layer_width_z-tol)) );

            bool on_edge_to_be_glued_to_sidewall=(
                    (potential_node_pt->x(1)<(y_min_first_sheet+tol))||
                    (potential_node_pt->x(1)>(y_max_first_sheet-tol))  );

            // Use set to avoid duplicates
            if (in_glued_strip||on_edge_to_be_glued_to_sidewall)
            {
                glue_node_set_pt.insert(potential_node_pt);
            }
        }

        // Transfer to vector
        for (std::set<Node*>::iterator it=glue_node_set_pt.begin();
             it!=glue_node_set_pt.end();it++)
        {
            glue_node_pt.push_back((*it));
        }

        // Glue specified nodes to co-located nodes in consituent mesh 3
        unsigned i_mesh_replace=3;
        Glued_mesh_pt->glue(glue_node_pt,i_mesh_replace);
    }





    // Second sheet glued to third sheet
    //----------------------------------
    {

        // Width of glued layer:
        double glue_layer_width_z=
                (z_max_first_sheet-z_min_first_sheet)/double(n_z_first_sheet);

        // We retain them on mesh 2 where the glued
        // nodes are on boundary 4
        unsigned i_mesh=2;
        unsigned b=4;
        unsigned n_node = Solid_mesh_pt[i_mesh]->nboundary_node(b);
        Vector<Node*> glue_node_pt;
        std::set<Node*> glue_node_set_pt;
        glue_node_pt.reserve(n_node);
        for(unsigned n=0;n<n_node;n++)
        {
            Node* potential_node_pt=Solid_mesh_pt[i_mesh]->boundary_node_pt(b,n);

            bool in_glued_strip=(
                    (potential_node_pt->x(2)>
                     (z_min_first_sheet+fold_height-glue_layer_width_z-tol)) );

            bool on_edge_to_be_glued_to_sidewall=(
                    (potential_node_pt->x(1)<(y_min_first_sheet+tol))||
                    (potential_node_pt->x(1)>(y_max_first_sheet-tol))  );

            // Use set to avoid duplicates
            if (in_glued_strip||on_edge_to_be_glued_to_sidewall)
            {
                glue_node_set_pt.insert(potential_node_pt);
            }
        }

        // Transfer to vector
        for (std::set<Node*>::iterator it=glue_node_set_pt.begin();
             it!=glue_node_set_pt.end();it++)
        {
            glue_node_pt.push_back((*it));
        }

        // Glue specified nodes to co-located nodes in consituent mesh 1
        unsigned i_mesh_replace=1;
        Glued_mesh_pt->glue(glue_node_pt,i_mesh_replace);

    }

    bool glue_sidewalls=true;
    if (glue_sidewalls)
    {
        // Left sidewall (at x=0)
        //-----------------------
        {
            bool on_left=true;
            std::set<SolidNode*> unique_nodes;
            get_unique_nodes_that_meet_sidewalls(on_left,unique_nodes);
            unsigned nnod=unique_nodes.size();
            Vector<Node*> glue_node_pt;
            glue_node_pt.reserve(nnod);
            std::set<SolidNode*>::iterator it;
            for (it=unique_nodes.begin();it!=unique_nodes.end();it++)
            {
                glue_node_pt.push_back((*it));
            }

            // Glue specified nodes to co-located nodes in consituent
            // mesh 4 (left sidewall)
            unsigned i_mesh_replace=4;
            Glued_mesh_pt->glue(glue_node_pt,i_mesh_replace);
        }

        // Right sidewall (at x=xmax)
        //---------------------------
        {
            bool on_left=false;
            std::set<SolidNode*> unique_nodes;
            get_unique_nodes_that_meet_sidewalls(on_left,unique_nodes);
            unsigned nnod=unique_nodes.size();
            Vector<Node*> glue_node_pt;
            glue_node_pt.reserve(nnod);
            std::set<SolidNode*>::iterator it;
            for (it=unique_nodes.begin();it!=unique_nodes.end();it++)
            {
                glue_node_pt.push_back((*it));
            }

            // Glue specified nodes to co-located nodes in
            // consituent mesh 5 (right sidewall)
            unsigned i_mesh_replace=5;
            Glued_mesh_pt->glue(glue_node_pt,i_mesh_replace);
        }
    }

    // Add to global mesh
    add_sub_mesh(Glued_mesh_pt);

    // Setup boundary conditions
    //--------------------------

    // hierher: use to apply other bcs; the minimally constrained one
    // isn't very good anyway (thought it does suppress rigid body modes)
    // because the entire resultant force from the pressure load has to be
    // balanced by a single point force which creates massive local deformations.
    bool minimal_constraint=true;


    // Suppress rigid body modes by pinning six displacements
    //--------------------------------------------------------
    // at lower-left, etc. nodes on
    //------------------------------
    // const y face at the glued end of the first sheet
    //-------------------------------------------------
    unsigned i_mesh=0;
    {
        unsigned b=4;
        SolidNode* ll_node_pt=Solid_mesh_pt[i_mesh]->boundary_node_pt(b,0);
        SolidNode* ur_node_pt=Solid_mesh_pt[i_mesh]->boundary_node_pt(b,0);
        SolidNode* ul_node_pt=Solid_mesh_pt[i_mesh]->boundary_node_pt(b,0);
        SolidNode* lr_node_pt=Solid_mesh_pt[i_mesh]->boundary_node_pt(b,0);

        unsigned n_node = Solid_mesh_pt[i_mesh]->nboundary_node(b);
        for(unsigned n=0;n<n_node;n++)
        {
            SolidNode* nod_pt=Solid_mesh_pt[i_mesh]->boundary_node_pt(b,n);
            double y=nod_pt->x(1);
            double z=nod_pt->x(2);
            if ((y<=ll_node_pt->x(1))&&(z<=ll_node_pt->x(2)))
            {
                ll_node_pt=nod_pt;
            }
            if ((y<=ul_node_pt->x(1))&&(z>=ul_node_pt->x(2)))
            {
                ul_node_pt=nod_pt;
            }
            if ((y>=ur_node_pt->x(1))&&(z>=ur_node_pt->x(2)))
            {
                ur_node_pt=nod_pt;
            }
            if ((y>=lr_node_pt->x(1))&&(z<=lr_node_pt->x(2)))
            {
                lr_node_pt=nod_pt;
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



    // Create traction meshes (attach to elements in still existing
    //-------------------------------------------------------------
    // original meshes)
    //-----------------

    // hierher: currently adding equal and opposite tractions on the exposed
    // (and glued together!) faces of the four sheets. Elements should be added
    // also the sidewalls and to the exposed top bit of the inner fold.
    // However, these will all be replaced by Twente's particle-based tractions
    // so they'll have to revisit this anyway.
    for (unsigned i_mesh=0;i_mesh<n_mesh;i_mesh++)
    {
        // Make new mesh
        Traction_mesh_pt[i_mesh] = new SolidMesh;

        // Create elements on x=1
        unsigned b=2;
        switch (i_mesh)
        {
            case 0:
            case 2:
                break;

            case 1:
            case 3:
                b=4;
                break;

            default:
                abort();
        }

        // Build the bastards
        create_traction_elements(Traction_mesh_pt[i_mesh],Solid_mesh_pt[i_mesh],b);

        // Add 'em
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
void SheetGlueProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{

    ofstream some_file;
    char filename[100];

    // Number of plot points
    unsigned npts = 5;

    // Output shape of deformed body
    sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    Solid_mesh_pt[0]->output(some_file,npts);
    Solid_mesh_pt[1]->output(some_file,npts);
    some_file.close();

    // Output traction
    //----------------
    sprintf(filename,"%s/traction%i.dat",doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    Traction_mesh_pt[0]->output(some_file,npts);
    Traction_mesh_pt[1]->output(some_file,npts);
    some_file.close();


    // Output glued mesh
    //------------------
    if (Glued_mesh_pt!=0)
    {
        sprintf(filename,"%s/glued_soln%i.dat",doc_info.directory().c_str(),
                doc_info.number());
        some_file.open(filename);
        Glued_mesh_pt->output(some_file,npts);
        some_file.close();
    }

}


//==================================================================
/// Run the problem
//==================================================================
template<class ELEMENT>
void SheetGlueProblem<ELEMENT>::run(const std::string &dirname)
{

    // Output
    DocInfo doc_info;

    // Set output directory
    doc_info.set_directory(dirname);

    // Step number
    doc_info.number()=0;

    // Doc initial configuration
    doc_solution(doc_info);
    doc_info.number()++;


    // Pressure increment
    double dp=0.0005;
    Global_Physical_Variables::P=0.0;

    //Parameter incrementation
    unsigned nstep=1000;
    if (CommandLineArgs::Argc!=1)
    {
        std::cout << "Validation -- only doing one step" << std::endl;
        nstep=1;
    }

    for(unsigned i=0;i<nstep;i++)
    {
        oomph_info << "Solving for p = "
                   << Global_Physical_Variables::P
                   << std::endl;

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
int main(int argc, char **argv)
{

    // Store command line arguments
    CommandLineArgs::setup(argc,argv);

    //Initialise physical parameters
    Global_Physical_Variables::E  = 2.1;
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
    SheetGlueProblem<QPVDElement<3,3> > problem;
    problem.run("./");
}





