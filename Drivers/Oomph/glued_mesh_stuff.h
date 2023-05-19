// hierher temporary header simply for inclusion of
// common code into multiple drivers; tidy up


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


//=========================================================================
/// Glued mesh; created by replacing pointers to spatially co-located nodes
/// of multiple constituent meshes. The constituent meshes stay alive
/// because they're used to identify mesh boundaries etc. Their boundary lookup
/// schemes etc. are updated to account for the "glued" nodes.
/// hierher currently only for SolidMesh; should be generalised before the
/// machinery is moved into src/generic.
/// hierher search for collocated nodes currently brute forced. Could be
/// optimised with cgal's tree-based search.
//=========================================================================
class GluedSolidMesh : public virtual SolidMesh {
    /// Vector of pointers to (still existing) constituent meshes;
    /// we retain them for access to their boundary enumerations
    Vector<SolidMesh *> Constituent_mesh_pt;

public:

    /// Constructor: Pass vector of source meshes. Co-lated nodes can be
    /// glued together (by re-allocating pointers)
    GluedSolidMesh(Vector<SolidMesh *> solid_mesh_pt) :
            Constituent_mesh_pt(solid_mesh_pt) {
        // Copy all nodes
        unsigned n_mesh = solid_mesh_pt.size();
        unsigned total_nnode = 0;
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++) {
            total_nnode += solid_mesh_pt[i_mesh]->nnode();
        }
        Node_pt.resize(total_nnode);

        unsigned node_count = 0;
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++) {
            unsigned nnode = solid_mesh_pt[i_mesh]->nnode();
            for (unsigned j = 0; j < nnode; j++) {
                Node_pt[node_count] = solid_mesh_pt[i_mesh]->node_pt(j);
                node_count++;
            }
        }

        // Copy all elements
        unsigned total_nel = 0;
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++) {
            total_nel += solid_mesh_pt[i_mesh]->nelement();
        }
        Element_pt.resize(total_nel);

        unsigned el_count = 0;
        for (unsigned i_mesh = 0; i_mesh < n_mesh; i_mesh++) {
            unsigned nel = solid_mesh_pt[i_mesh]->nelement();
            for (unsigned e = 0; e < nel; e++) {
                Element_pt[el_count] = solid_mesh_pt[i_mesh]->element_pt(e);
                el_count++;
            }
        }
    }


    /// Glue specified nodes to co-located nodes
    /// in consituent mesh i_mesh_replace. i.e pointers to
    /// nodes in constituent mesh i_mesh_replace are replaced by pointers
    /// to colocated nodes (in glue_node_pt). Original nodes are then deleted
    /// and the boundary lookup schemes of constituent mesh i_mesh_replace
    /// are updated.
    void glue(const Vector<Node *> glue_node_pt,
              const unsigned &i_mesh_replace) {
        double t_start = TimingHelpers::timer();

        // Make input parameter
        double tol = 1.0e-10;

        // Map from original to replacement node:
        // Usage: replacement_node_pt[replaced_node_pt]=nod_pt
        std::map<Node *, Node *> replacement_node_pt;

        // Brute force it:
        //================
        unsigned nnod_glue = glue_node_pt.size();
        for (unsigned j = 0; j < nnod_glue; j++) {
            Node *nod_pt = glue_node_pt[j];
            unsigned dim = nod_pt->ndim();
            Node *node_to_be_replaced_pt = 0;
            unsigned jj_replace_index = 0;
            unsigned nnod_candidate = Constituent_mesh_pt[i_mesh_replace]->nnode();
            for (unsigned jj = 0; jj < nnod_candidate; jj++) {
                Node *nod_alt_pt = Constituent_mesh_pt[i_mesh_replace]->node_pt(jj);
                double dist_squared = 0.0;
                for (unsigned i = 0; i < 3; i++) {
                    dist_squared +=
                            (nod_alt_pt->x(i) - nod_pt->x(i)) *
                            (nod_alt_pt->x(i) - nod_pt->x(i));
                }
                if (sqrt(dist_squared) < tol) {
                    node_to_be_replaced_pt = nod_alt_pt;
                    jj_replace_index = jj;
                    break;
                }
            }
            if (node_to_be_replaced_pt == 0) {
                oomph_info << "ERROR: Not found a replacement node for node at ";
                for (unsigned i = 0; i < dim; i++) {
                    oomph_info << nod_pt->x(i) << " ";
                }
                oomph_info << std::endl;
                // hierher throw proper error
                abort();
            }

            // Replace node in constituent mesh
            Node *replaced_node_pt =
                    Constituent_mesh_pt[i_mesh_replace]->node_pt(jj_replace_index);
            dynamic_cast<Mesh *>(Constituent_mesh_pt[i_mesh_replace])->
                    node_pt(jj_replace_index) = nod_pt;
            replacement_node_pt[replaced_node_pt] = nod_pt;

            // Replace node in elements
            unsigned nel = Constituent_mesh_pt[i_mesh_replace]->nelement();
            for (unsigned e = 0; e < nel; e++) {
                FiniteElement *fe_pt =
                        Constituent_mesh_pt[i_mesh_replace]->finite_element_pt(e);
                unsigned nod_el = fe_pt->nnode();
                for (unsigned j_in_el = 0; j_in_el < nod_el; j_in_el++) {
                    if (fe_pt->node_pt(j_in_el) == node_to_be_replaced_pt) {
                        fe_pt->node_pt(j_in_el) = nod_pt;
                    }
                }
            }

            // Replace node in original mesh's boundary lookup scheme
            unsigned nb = Constituent_mesh_pt[i_mesh_replace]->nboundary();
            for (unsigned b = 0; b < nb; b++) {
                unsigned nnod = Constituent_mesh_pt[i_mesh_replace]->nboundary_node(b);
                for (unsigned j = 0; j < nnod; j++) {
                    Node *potentially_replaced_node_pt =
                            Constituent_mesh_pt[i_mesh_replace]->boundary_node_pt(b, j);
                    std::map<Node *, Node *>::iterator it;
                    it = replacement_node_pt.find(potentially_replaced_node_pt);
                    if (it == replacement_node_pt.end()) {
                        //oomph_info << "Node in boundary lookup scheme was not replaced\n";
                    } else {
                        // oomph_info << "Replace node "
                        //            << potentially_replaced_node_pt << " at: "
                        //            << (it->second)->x(0) << " "
                        //            << (it->second)->x(1) << " "
                        //            << (it->second)->x(2) << " "
                        //            << "with " << it->second
                        //            << std::endl;

                        // hierher: SolidMesh needs an overload of the read/write version
                        // to this function
                        dynamic_cast<Mesh *>(Constituent_mesh_pt[i_mesh_replace])->
                                boundary_node_pt(b, j) = it->second;
                    }
                }
            }


            // Remove it from the glued mesh's compound Node_pt vector:
            Vector<Node *> tmp_node_pt(Node_pt);
            unsigned nnod = tmp_node_pt.size();
            Node_pt.resize(nnod - 1);
            unsigned count = 0;
            for (unsigned jj = 0; jj < nnod; jj++) {
                if (tmp_node_pt[jj] != node_to_be_replaced_pt) {
                    Node_pt[count] = tmp_node_pt[jj];
                    count++;
                }
            }

            // Now kill node
            // hierher REINSTATE deletion
            delete node_to_be_replaced_pt;
        }

        // oomph_info << "After gluing glued mesh has " << Node_pt.size()
        //            << " nodes" << std::endl;


        // Tell us how bad it was; do you fancy implementing the optimised
        // search scheme?
        double t_end = TimingHelpers::timer();
        oomph_info << "Time for (inefficient!) mesh gluing: "
                   << t_end - t_start << std::endl;


    }
};



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


//=========================================================================
/// Simple cubic mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class SolidCubicMesh : public virtual SimpleCubicMesh<ELEMENT>,
                       public virtual SolidMesh {

public:

    /// Constructor: Specify half widths
    SolidCubicMesh(const unsigned &nx, const unsigned &ny,
                   const unsigned &nz,
                   const double &a,
                   const double &b,
                   const double &c,
                   TimeStepper *time_stepper_pt =
                   &Mesh::Default_TimeStepper) :
            SimpleCubicMesh<ELEMENT>(nx, ny, nz, -a, a, -b, b, -c, c, time_stepper_pt),
            SolidMesh() {
        //Assign the initial lagrangian coordinates
        set_lagrangian_nodal_coordinates();
    }


    /// Constructor:
    SolidCubicMesh(const unsigned &nx,
                   const unsigned &ny,
                   const unsigned &nz,
                   const double &x_min,
                   const double &x_max,
                   const double &y_min,
                   const double &y_max,
                   const double &z_min,
                   const double &z_max,
                   TimeStepper *time_stepper_pt =
                   &Mesh::Default_TimeStepper) :
            SimpleCubicMesh<ELEMENT>(nx, ny, nz,
                                     x_min, x_max,
                                     y_min, y_max,
                                     z_min, z_max,
                                     time_stepper_pt),
            SolidMesh() {
        //Assign the initial lagrangian coordinates
        set_lagrangian_nodal_coordinates();
    }

    /// Empty Destructor
    virtual ~SolidCubicMesh() {}

};


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
