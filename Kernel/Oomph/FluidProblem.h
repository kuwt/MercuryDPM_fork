//
// Created by mitchel on 11/7/22.
//

#ifndef MERCURYDPM_FLUIDPROBLEM_H
#define MERCURYDPM_FLUIDPROBLEM_H

// Generic oomph-lib headers
#include "generic.h"
#include "Math/Vector.h"
#include <array>
#include <cstdio>
#include "OomphHelpers.h"
#include "mesh.h"
#include "meshes/simple_cubic_mesh.h"
#include "navier_stokes/navier_stokes.h"

using namespace oomph;

template<class ELEMENT_TYPE>
class FluidProblem : public Problem
{
protected:
    typedef ELEMENT_TYPE ELEMENT;
    
    RefineableMeshBase* Fluid_mesh_pt = nullptr;
    double Re_ = 0.0;
    double ReSt_ = 0.0;
    double Re_InvFr_ = 0.0;
    
    // Fluid properties
    double fluidDensity = 0.0;
    double fluidKinematicViscosity = 0.0;
    double fluidDynamicViscosity = 0.0;
    int adaptEveryNFluidTimesteps_ = 0;
    
public:
    FluidProblem()
    {
        logger(INFO, "FluidProblem constructor, setting timestepper to BDF<2>");
        
        // Set the timestepper
        add_time_stepper_pt(new BDF<2>);
        
        Max_newton_iterations = 200;
        Newton_solver_tolerance = 1e-10;
        Max_residuals = constants::inf;
    }
    
    void refineMesh()
    {
        if (adaptEveryNFluidTimesteps_ <= 0)
        {
            logger(WARN, "adaptEveryNFluidTimesteps_ is not set to value > 0, mesh will never be refined even though element is refineable");
        }
        else
        {
            if (static_cast<int>(this->time_pt()->time() / getOomphTimeStep()) % getAdaptEveryNFluidTimesteps() == 0)
            {
                logger(DEBUG, "Call to adapt()");
                adapt();
            }
        }
    }
    
    class RefineableFluidCubicMesh : public virtual RefineableSimpleCubicMesh<ELEMENT>, public virtual RefineableMeshBase
    {
    public:
        /// \short Constructor:
        // nx, ny, nz: number of elements in the x, y, and z directions
        // xMax, yMax, zMax: dimensions of the cube (assume the center of the cube is at the origin)
        // timeStepper: defaults to Steady.
        RefineableFluidCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                       const double& xMin, const double& xMax, const double& yMin,
                       const double& yMax, const double& zMin, const double& zMax,
                       TimeStepper* time_stepper_pt) :
                SimpleCubicMesh<ELEMENT>(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt),
                RefineableSimpleCubicMesh<ELEMENT>(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt),
                RefineableMeshBase()
        {
        
        }
    };
    
    // Create a solid cubic mesh
    void setRefineableFluidCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                           const double& xMin, const double& xMax, const double& yMin,
                           const double& yMax, const double& zMin, const double& zMax)
    {
        // Create mesh
        fluid_mesh_pt() = new RefineableFluidCubicMesh(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt());
        
        Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
        fluid_mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
        
        logger(INFO, "Created %x%x% refineable cubic mesh on domain [%,%]x[%,%]x[%,%]",
               nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax);
    }
    
    class FluidCubicMesh : public virtual SimpleCubicMesh<ELEMENT>
    {
    public:
        /// \short Constructor:
        // nx, ny, nz: number of elements in the x, y, and z directions
        // xMax, yMax, zMax: dimensions of the cube (assume the center of the cube is at the origin)
        // timeStepper: defaults to Steady.
        FluidCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                                 const double& xMin, const double& xMax, const double& yMin,
                                 const double& yMax, const double& zMin, const double& zMax,
                                 TimeStepper* time_stepper_pt) :
                SimpleCubicMesh<ELEMENT>(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt)
        {
        
        }
    };
    
    // Create a solid cubic mesh
    void setFluidCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                                     const double& xMin, const double& xMax, const double& yMin,
                                     const double& yMax, const double& zMin, const double& zMax)
    {
        // Create mesh
        fluid_mesh_pt() = new FluidCubicMesh(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt());
        
        Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
        fluid_mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
        
        logger(INFO, "Created %x%x% cubic mesh on domain [%,%]x[%,%]x[%,%]",
               nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax);
    }
    
    virtual void actionsBeforeOomphTimeStep() {}
    
    /// Get function for the fluid mesh pointer
    RefineableMeshBase*& fluid_mesh_pt() { return Fluid_mesh_pt; }
    /// Get function for the solid mesh pointer
    RefineableMeshBase* const & fluid_mesh_pt() const { return Fluid_mesh_pt; }
    
    
    // Functions regarding scaling
    void setReynoldsNumber(const double& Re) {Re_ = Re;}
    double getReynoldsNumber() {return Re_;}
    void setReynoldsStrouhalNumber(const double& ReSt) {ReSt_ = ReSt;}
    double getReynoldsStrouhalNumber() {return ReSt_;}
    void setReynolds_InverseFroudeNumber(const double& Re_InvFr) {Re_InvFr_ = Re_InvFr;}
    double getReynolds_InverseFroudeNumber() {return Re_InvFr_;}
    
    void setAdaptEveryNFluidTimesteps(const int nFluidTimesteps) {adaptEveryNFluidTimesteps_ = nFluidTimesteps;}
    int getAdaptEveryNFluidTimesteps() {return adaptEveryNFluidTimesteps_;}
    
    void pinBC(unsigned iBoundary_)
    {
        logger(INFO,"Pinning boundary %",iBoundary_);
        
        unsigned long int nNode = fluid_mesh_pt()->nboundary_node(iBoundary_);
    
        for (unsigned iNode = 0; iNode < nNode; iNode++)
        {
            //Pin u and v and w
            fluid_mesh_pt()->boundary_node_pt(iBoundary_, iNode)->pin(0);
            fluid_mesh_pt()->boundary_node_pt(iBoundary_, iNode)->pin(1);
            fluid_mesh_pt()->boundary_node_pt(iBoundary_, iNode)->pin(2);
        }
    };
    
    void setBC(unsigned iBoundary_, unsigned dof_, double value_)
    {
        logger(INFO,"Setting boundary %, dof % at value %",iBoundary_,dof_,value_);
    
        unsigned long int nNode = fluid_mesh_pt()->nboundary_node(iBoundary_);
        for (unsigned iNode = 0; iNode < nNode; iNode++)
        {
            //logger(INFO,"Value is not pinned == %",!(fluid_mesh_pt()->boundary_node_pt(iBoundary_,iNode)->is_pinned(dof_)));
            logger.assert_always(fluid_mesh_pt()->boundary_node_pt(iBoundary_,iNode)->is_pinned(dof_),"The dof you are tryin to set is not pinned");
            fluid_mesh_pt()->boundary_node_pt(iBoundary_, iNode)->set_value(dof_,value_);
        }
    }
    
    /// set function for name_
    void setName(const std::string& name)
    {
        name_ = name;
        logger(INFO, "Name: %", name_);
    }
    
    /// set function for time step
    void setOomphTimeStep(double dt)
    {
        time_pt()->initialise_dt(dt);
    }
    double getOomphTimeStep()
    {
        return time_pt()->dt();
    }
    
    void prepareForSolve()
    {
        logger.assert_always(fluidDensity>0,"fluidDensity not initialised");
    
        logger.assert_always(fluid_mesh_pt(), "Set fluid mesh via e.g. setRefineableFluidCubicMesh(..)");
    
        logger(INFO, "Assign re, re_st, re_inv_fr");
        for (unsigned i = 0; i < fluid_mesh_pt()->nelement(); i++)
        {
            //Cast to a fluid element
            ELEMENT* el_pt = dynamic_cast<ELEMENT*>(fluid_mesh_pt()->element_pt(i));
            el_pt->re_pt()= &Re_;
            el_pt->re_st_pt() = &ReSt_;
            el_pt->re_invfr_pt() = &Re_InvFr_;
        }
    
        add_sub_mesh(fluid_mesh_pt());
        build_global_mesh();
        
        // Pin the redundant fluid pressures (if any)
        ELEMENT::pin_redundant_nodal_pressures(fluid_mesh_pt()->element_pt());
        logger(INFO, "Pinned redundant nodal fluid pressures");
    
        // Attach the boundary conditions to the mesh
        unsigned n_eq = assign_eqn_numbers();
        logger(INFO, "Assigned % equation numbers", n_eq);
    }
    
    /// Doc the solution
    void doc_solution(DocInfo& doc_info);
    /// Doc the solution paraview format
    void doc_paraview(DocInfo& doc_info);
    /// Doc the voidage
    void doc_voidage(DocInfo& doc_info);
    /// Doc the elements
    void doc_element(DocInfo& doc_info);
    
    // Setters and getters for fluid-properties
    void setFluidDensity(double fluidDensity_) {fluidDensity = fluidDensity_;}
    double getFluidDensity() {return fluidDensity;}
    void setFluidDynamicViscosity(double fluidDynamicViscosity_) {fluidDynamicViscosity = fluidDynamicViscosity_;}
    double getFluidDynamicViscosity() {return fluidDynamicViscosity;}
    void setFluidKinematicViscosity(double fluidKinematicViscosity_) {fluidKinematicViscosity = fluidKinematicViscosity_;}
    double getFluidKinematicViscosity() {return fluidKinematicViscosity;}

protected:
    
    std::string name_;
    DocInfo doc_info;
    
};

//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void FluidProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
    if (doc_info.directory().empty())
    {
        logger(ERROR,"Directory for oomph-output is not set, please use doc_info.set_directory(\"path\")");
    }
    
    std::ofstream some_file;
    char filename[100];
    
    // Number of plot points
    unsigned npts;
    npts=5;
    
    // Output solution
    sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),doc_info.number());
    some_file.open(filename);
    fluid_mesh_pt()->output(some_file,npts);
    some_file.close();
} //end doc_solution


//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void FluidProblem<ELEMENT>::doc_paraview(DocInfo& doc_info)
{
    if (doc_info.directory().empty())
    {
        logger(ERROR,"Directory for oomph-output is not set, please use doc_info.set_directory(\"path\")");
    }
    
    std::ofstream some_file;
    char filename[100];
    
    // Number of plot points
    unsigned npts;
    npts=5;
    
    // Output solution
    sprintf(filename,"%s/soln%i.pvd",doc_info.directory().c_str(),doc_info.number());
    some_file.open(filename);
    fluid_mesh_pt()->output_paraview(some_file,npts);
    some_file.close();
} //end doc_solution

//==============start_doc_void======================================
/// Doc the voidage ///FIXME Needs function pointer to output_voidage_byEl(some_file, npts) which only exists for AJ equations
//==================================================================
template<class ELEMENT>
void FluidProblem<ELEMENT>::doc_voidage(DocInfo& doc_info)
{
    if (doc_info.directory().empty())
    {
        logger(ERROR,"Directory for oomph-output is not set, please use doc_info.set_directory(\"path\")");
    }

    std::ofstream some_file;
    char filename[100];
    
    // number of plot points
    unsigned npts = 1;
    
    // Output solution if using get_voidage_byEl
    sprintf(filename,"%s/voidagen%i.dat",doc_info.directory().c_str(),doc_info.number());
    some_file.open(filename);
    //fluid_mesh_pt()->output_voidage_byEl(some_file,npts); //FIXME output_voidage_byEl function pointer not set for non-AJ elements
    some_file.close();
} //end doc_voidage

//==============start_doc_element======================================
/// Doc the element data
//==================================================================
template<class ELEMENT>
void FluidProblem<ELEMENT>::doc_element(DocInfo& doc_info)
{
    if (doc_info.directory().empty())
    {
        logger(ERROR,"Directory for oomph-output is not set, please use doc_info.set_directory(\"path\")");
    }
    
    std::ofstream some_file;
    char filename[100];
    
    // number of plot points
    unsigned npts = 1;
    
    // Output solution
    sprintf(filename,"%s/elements%i.dat",doc_info.directory().c_str(),doc_info.number());
    some_file.open(filename);
    fluid_mesh_pt()->output(some_file,npts);
    some_file.close();
} //end doc_voidage

#endif //MERCURYDPM_FLUIDPROBLEM_H
