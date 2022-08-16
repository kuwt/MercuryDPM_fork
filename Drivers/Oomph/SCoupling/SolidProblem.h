//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://MercuryDPM.org/Team>.
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
 * For the material model, see
 * http://oomph-lib.maths.man.ac.uk/doc/solid/solid_theory/html/index.html
 **/

#ifndef SOLID_PROBLEM_H
#define SOLID_PROBLEM_H

// Generic oomph-lib headers
#include "generic.h"
#include "solid.h"
#include "constitutive.h"
#include "Math/Vector.h"
#include <array>
#include <cstdio>
#include "Coupling/AnisotropicHookean.h"
#include "Coupling/RefineableQDPVDElement.h"

//The mesh
#include "meshes/simple_cubic_mesh.h"
#include "meshes/tetgen_mesh.h"

// The element types needed for surface and volume coupling
#include "Coupling/SurfaceCoupledElement.h"
#include "Coupling/VolumeCoupledElement.h"

using namespace oomph;

/**
 * Base class for (surface-coupled) problems with a solid body.
 *
 * Assumptions:
 *  - The element type is derived from QPVDElement, i.e. a solid.
 *  - The mesh type is derived from SolidMesh.
 *  - The constitutive law is derived from Hookean, with parameters elasticModulus, poissonRatio, density.
 *
 * The solid properties (elasticModulus, poissonRatio, density, constitutive_law_pt, body_force_fct, etc.) are stored in member variables, not in globals as typical in oomph-lib, and accessed by setters and getters.
 *
 * Additional functionality:
 * - setIsPinned allows one to define the pinned nodes using a function.
 */
template<class ELEMENT>
class SolidProblem : public Problem
{
protected:
    //define element type (should be done by template)
    //typedef QPVDElement<3, 2> ELEMENT;

    // name of output files (user-defined)
    std::string name_;

    /// Elastic modulus (set via setSolid)
    double elasticModulus_ = 0;

    /// Poisson's ratio (set via setSolid)
    double poissonRatio_ = 0;

    /// Density
    double density_ = 0;

    /// Pointer to the body force function
    void (* body_force_fct)(const double& time, const Vector<double>& xi, Vector<double>& b) = nullptr;

    /// Pointer to constitutive law (should be set in constructor)
    ConstitutiveLaw* constitutive_law_pt  = nullptr;

    /// Pointer to solid mesh
    SolidMesh* Solid_mesh_pt = nullptr;

    /// Pointer to mesh of traction elements
    SolidMesh* Traction_mesh_pt = nullptr;

    /// Function to determine which nodal positions are pinned
    std::function<bool(SolidNode*, unsigned)> isPinned_;

public:

    /// Empty constructor:
    SolidProblem()
    {
        // Set Newton solver tolerance and maximum number of Newton iterations
        Max_newton_iterations = constants::unsignedMax;
        Newton_solver_tolerance = 1e-10;
        Max_residuals = constants::inf;

        // Set constitutive law
        constitutive_law_pt = new GeneralisedHookean(&poissonRatio_, &elasticModulus_);

        // Allocate the timestepper
        add_time_stepper_pt(new Newmark<2>);

        //build_mesh();

        // setup boundary conditions
        //pin_and_assign_eqn_numbers();
    };

    /// set function for name_
    void setName(const std::string& name)
    {
        name_ = name;
        logger(INFO, "Name: %", name_);
    }

    /// set function for elasticModulus_
    void setElasticModulus(double elasticModulus)
    {
        elasticModulus_ = elasticModulus;
        logger(INFO, "Elastic Modulus: %", elasticModulus_);
    }

    /// set function for poissonRatio_
    void setPoissonRatio(double poissonRatio)
    {
        poissonRatio_ = poissonRatio;
        logger(INFO, "Poisson Ratio: %", poissonRatio_);
    }

    /// set function for density_
    void setDensity(double density)
    {
        density_ = density;
        logger(INFO, "Density: %", density_);
    }

    /// set function for body_force_pt
    void setBodyForceAsGravity()
    {
        // define a static body force
        static double& Density = density_;
        body_force_fct = [](const double& time, const Vector<double>& xi, Vector<double>& b) {
            b[0] = 0.0;
            b[1] = 0.0;
            b[2] = -9.8 * Density;
        };
    }

    /// set function for isPinned_
    void setIsPinned(std::function<bool(SolidNode*, unsigned)> isPinned)
    {
        isPinned_ = std::move(isPinned);
        logger(INFO, "Setting which positions on which nodes are pinned");
    }

    // set is_pinned such that a certain boundary is pinned
    void pinBoundary(unsigned b) {
        isPinned_ = [b](SolidNode *n, unsigned d) {
            return n->is_on_boundary(b);
        };
    }

    /// Sets the dissipation coefficient for all elements.
    std::enable_if<std::is_base_of<RefineableQDPVDElement<3,2>, ELEMENT>::value, void> setDissipation(double dissipation)
    {
        logger(INFO,"Adding dissipation %",dissipation);
        for (int i = 0; i < solid_mesh_pt()->nelement(); ++i)
        {
            dynamic_cast<RefineableQDPVDElement<3,2>*>(solid_mesh_pt()->element_pt(i))->setDissipation(dissipation);
        }
    }

    /// set function for Newton_solver_tolerance
    void setNewtonSolverTolerance(double Newton_solver_tolerance)
    {
        this->Newton_solver_tolerance = Newton_solver_tolerance;
    }

    /// set function for Max_newton_iterations
    void setMaxNewtonIterations(unsigned Max_newton_iterations)
    {
        this->Max_newton_iterations = Max_newton_iterations;
    }

    /// set function for time step
    void setOomphTimeStep(double dt) {
        time_pt()->initialise_dt(dt);
    }

    /// get function for time step
    double getOomphTimeStep () const {
        return time_pt()->dt();
    }

    /// get function for current time
    double getOomphTime () const {
        return time_pt()->time();
    }

    /// Get function for the solid mesh pointer
    SolidMesh*& solid_mesh_pt() { return Solid_mesh_pt; }

    /// Get function for the traction mesh pointer
    SolidMesh*& traction_mesh_pt() { return Traction_mesh_pt; }

    /// Computes the domain size (min/max of the nodal positions in x/y/z)
    void getDomainSize(std::array<double, 3>& min, std::array<double, 3>& max)
    {
        min[0] = min[1] = min[2] = constants::inf;
        max[0] = max[1] = max[2] = -constants::inf;
        logger.assert_always(solid_mesh_pt(),"mesh not found");
        for (unsigned i = 0; i < solid_mesh_pt()->nnode(); i++)
        {
            const auto n = solid_mesh_pt()->node_pt(i);
            for (int j = 0; j < 3; ++j)
            {
                min[j] = std::min(min[j], n->xi(j));
                max[j] = std::max(max[j], n->xi(j));
            }
        }
    }

    /**
     * - Asserts all parameters are set
     * - Assign pointers for all elements
     * - Builds global mesh
     * - Pins nodes
     * - Pins solid pressure
     * - Assigns equation numbers
     *
     * This function should be called before solve(), after creating the mesh, defining body force and constitutive law.
     */
    void prepareForSolve()
    {
        // check certain values are set
        logger.assert_always(!name_.empty(), "Set name via setName(..)");
        logger.assert_always(elasticModulus_>0, "Set elastic modulus via setElasticModulus(..)");
        logger.assert_always(density_>0, "Set density via setDensity(..)");
        logger.assert_always(solid_mesh_pt(), "Set solid mesh via e.g. setSolidCubicMesh(..)");

        // Assign constitutive_law_pt and body_force_fct_pt of each element
        for (unsigned i = 0; i < solid_mesh_pt()->nelement(); i++)
        {
            //Cast to a solid element
            ELEMENT* el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(i));
            // Set the constitutive law
            el_pt->constitutive_law_pt() = constitutive_law_pt;
            // Set body force
            el_pt->body_force_fct_pt() = body_force_fct;
            // Set density
            el_pt->lambda_sq_pt() = &density_;
        }

        // Construct the traction element mesh
        // Traction_mesh_pt = new SolidMesh;
        // create_traction_elements();

        //Build combined "global" mesh
        add_sub_mesh(solid_mesh_pt());
        if (traction_mesh_pt()) {
            add_sub_mesh(traction_mesh_pt());
            logger(INFO, "Built global mesh from solid mesh and traction mesh");
        } else {
            logger(INFO, "Built global mesh from solid mesh");
        }
        build_global_mesh();

        //Pin nodes
        if (isPinned_)
        {
            for ( unsigned n = 0; n < solid_mesh_pt()->nnode(); n++ )
            {
                SolidNode* n_pt = solid_mesh_pt()->node_pt( n );
                //Pin all nodes
                for ( unsigned i = 0; i < 3; i++ )
                {
                    if ( isPinned_( n_pt, i ) )
                    {
                        n_pt->pin_position( i );
                    }
                    else
                    {
                        n_pt->unpin_position( i );
                    }
                }
            }
        }
        countPinned();

        // Pin the redundant solid pressures (if any)
        PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
            solid_mesh_pt()->element_pt());
        logger(INFO, "Pinned redundant nodal solid pressures");

        // Attach the boundary conditions to the mesh
        unsigned n_eq = assign_eqn_numbers();
        logger(INFO, "Assigned % equation numbers", n_eq);
    }

    /// returns statistics about pinned nodes to the console
    void countPinned()
    {
        // count pinned
        std::array<unsigned,3> countPinned {0,0,0};
        unsigned countAll = 0;
        for (unsigned n = 0; n < solid_mesh_pt()->nnode(); n++)
        {
            for (unsigned i = 0; i < 3; i++)
            {
                countPinned[i] += solid_mesh_pt()->node_pt(n)->position_is_pinned(i);
                countAll++;
            }
        }
        unsigned countPinnedAll = countPinned[0] + countPinned[1] +countPinned[2];
        logger(INFO, "Pinned % of % positions (% free): % in x, % in y, % in z", countPinnedAll, countAll, countAll - countPinnedAll, countPinned[0], countPinned[1], countPinned[2]);
    }

    virtual void actionsBeforeOomphTimeStep() {}

    /**
     * Solves an unsteady problem.
     */
    void solveUnsteady(double timeMax, double dt, unsigned saveCount = 10)
    {
        std::cout << "Solving oomph with dt=" << dt << " until timeMax=" << timeMax << std::endl;

        // Setup initial conditions. Default is a non-impulsive start
        assign_initial_values_impulsive(dt);

        // This is the main loop over advancing time
        double& time = time_stepper_pt()->time();
        unsigned count = 0;
        const unsigned countMax = std::ceil(timeMax/dt);
        while (time < timeMax)
        {
            logger(INFO, "Time %s of %s (% of %)",
                    //logger(INFO, "\n\033[1;33mTime %s of %s\033[0m\n",
                    time * Global_Physical_Variables::timeScale,
                    timeMax * Global_Physical_Variables::timeScale, count, countMax);
            actionsBeforeOomphTimeStep();
            // solve the oomphProb for one time step (this also increments time)
            unsteady_newton_solve(dt);
            // increase count
            count++;
            // write outputs of the oomphProb
            if (count % saveCount == 0 or time + dt > timeMax) {
                writeToVTK();
            }
            // abort if problem
            if (Max_res.back()==0) {
                logger(WARN,"Maximum residual is 0; aborting time loop");
                break;
            }
        }
        saveSolidMesh();
    }

    /**
     * Solves a steady problem.
     */
    void solveSteady()
    {
        newton_solve();
        writeToVTK();
        saveSolidMesh();
    }

    /**
     * Removes old output files (vtu) with the same name as this problem.
     */
    void removeOldFiles() const
    {
        for (int i = 0; true; ++i)
        {
            std::string fileName = name_ + "FEM_" + std::to_string(i) + ".vtu";
            if (remove(fileName.c_str())) break;
            std::cout << "Deleted " << fileName << '\n';
        }
    }

    /**
     * Returns the x value at a certain xi value
     */
    void get_x(const Vector<double>& xi, Vector<double>& x)
    {
    #ifdef OOMPH_HAS_MPI
        logger(INFO,"get_x does not work with MPI");
    #else
        Vector<double> s(3);
        GeomObject* geom_obj_pt = nullptr;
        const unsigned long nelement = solid_mesh_pt()->nelement();
        for (unsigned long i = 0; i < nelement; i++)
        {
            auto el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(i));
            el_pt->locate_zeta(xi, geom_obj_pt, s);
            if (geom_obj_pt)
            {
                //logger(INFO,"Point % % % is in element % at % % %",
                //       xi[0],xi[1],xi[2],i,s[0], s[1], s[2]);
                el_pt->interpolated_x(s, x); //deformed coordinate
                return;
            }
        }
        logger(ERROR, "x(xi) could not be found");
    #endif
    }

    /**
     * Returns the difference between the x and xi value in dimension d at a certain xi value
     */
    double getDeflection(Vector<double> xi, unsigned d)
    {
        Vector<double> x(3);
        get_x(xi, x);
        return x[d] - xi[d];
    }

    // Boundary types (from cubic mesh)
    enum Boundary : unsigned
    {
        Z_MIN = 0, Y_MIN = 1, X_MAX = 2, Y_MAX = 3, X_MIN = 4, Z_MAX = 5
    };

    // Simple cubic mesh upgraded to become a solid mesh
    class SolidCubicMesh : public virtual RefineableSimpleCubicMesh<ELEMENT>, public virtual SolidMesh
    {
    public:
        /// \short Constructor:
        // nx, ny, nz: number of elements in the x, y, and z directions
        // xMax, yMax, zMax: dimensions of the cube (assume the center of the cube is at the origin)
        // timeStepper: defaults to Steady.
        SolidCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                       const double& xMin, const double& xMax, const double& yMin,
                       const double& yMax, const double& zMin, const double& zMax,
                       TimeStepper* time_stepper_pt) :
                SimpleCubicMesh<ELEMENT>(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt),
                RefineableSimpleCubicMesh<ELEMENT>(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt),
                SolidMesh()
        {
            //Assign the initial lagrangian coordinates
            set_lagrangian_nodal_coordinates();
        }
    };

    // Create a solid cubic mesh
    void setSolidCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                           const double& xMin, const double& xMax, const double& yMin,
                           const double& yMax, const double& zMin, const double& zMax)
    {
        // Create mesh
        solid_mesh_pt() = new SolidCubicMesh(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt());

        logger(INFO, "Created %x%x% cubic mesh on domain [%,%]x[%,%]x[%,%]",
               nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax);
    }

    void saveSolidMesh()
    {
        auto solid_cubic_mesh_pt = dynamic_cast<SolidCubicMesh*>(solid_mesh_pt());
        if (not solid_cubic_mesh_pt) {
            logger(INFO,"Solid mesh can only be stored if based on cubic constructor");
            return;
        }

        std::ofstream mesh(name_ + ".mesh");
        mesh << solid_cubic_mesh_pt->nx() << ' ';
        mesh << solid_cubic_mesh_pt->ny() << ' ';
        mesh << solid_cubic_mesh_pt->nz() << '\n';
        for (int i = 0; i < solid_mesh_pt()->nnode(); ++i)
        {
            SolidNode* n = solid_mesh_pt()->node_pt(i);
            for (int j = 0; j < 3; ++j)
            {
                mesh << n->xi(j) << ' ' << n->x(j) << ' ' << n->position_is_pinned(j) << ' ';
            }
            mesh << '\n';
        }
        logger(INFO, "Saved %x%x% mesh to %.mesh",solid_cubic_mesh_pt->nx(), solid_cubic_mesh_pt->ny(), solid_cubic_mesh_pt->nz(),name_);
    }

    void loadSolidMesh(std::string name)
    {
        //if (not std::is_base_of<QElement<3, 2>, ELEMENT>::value) {
        //    logger(ERROR, "This kind of mesh requires a QElement.");
        //}

        std::ifstream mesh(name);
        logger.assert_always(mesh.good(),"Mesh file % could not be opened",name);
        unsigned nx, ny, nz;
        mesh >> nx >> ny >> nz;
        logger(INFO, "Loaded %x%x% cubic mesh from %", nx, ny, nz, name);
        solid_mesh_pt() = new SolidCubicMesh(nx, ny, nz, 0, 1, 0, 1, 0, 1, time_stepper_pt());

        double xi, x;
        bool pin;
        for (int i = 0; i < solid_mesh_pt()->nnode(); ++i)
        {
            SolidNode* n = solid_mesh_pt()->node_pt(i);
            for (int j = 0; j < 3; ++j)
            {
                mesh >> xi >> x >> pin;
                n->xi(j) = xi;
                n->x(j) = x;
                if (pin) n->pin_position(j);
            }
        }
    }

    // stores results in vtk file
    void writeToVTK()
    {
        //set local coordinates list
        std::vector<std::vector<double>> sList0;
        // order of nodes
        std::vector<unsigned> nList;
        // vtkFormat for given shape
        // https://kitware.github.io/vtk-examples/site/VTKFileFormats/
        unsigned vtkFormat;

        // define differently for tetrahedral and hexahedral (quad) elements
        if (std::is_base_of<SolidTElement<3, 2>, ELEMENT>::value)
        {
            vtkFormat = 10; // TETRA

            sList0 = {
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1},
                {0, 0, 0}
            };

            nList = {0, 1, 2, 3};
        }
        else if (std::is_base_of<SolidQElement<3, 2>, ELEMENT>::value)
        {
            vtkFormat = 12; // HEXAHEDRON

            sList0 = {
                {-1, -1, -1},
                {-1, -1, +1},
                {-1, +1, +1},
                {-1, +1, -1},
                {+1, -1, -1},
                {+1, -1, +1},
                {+1, +1, +1},
                {+1, +1, -1}
            };
            // order of nodes
            nList = {0, 4, 6, 2, 1, 5, 7, 3};
        } else {
            logger(ERROR,"Element type unknown");
        }

        // convert to Vector
        std::vector<Vector<double>> sList;
        for (auto s0 : sList0)
        {
            Vector<double> s(3);
            s[0] = s0[0];
            s[1] = s0[1];
            s[2] = s0[2];
            sList.push_back(s);
        }

        // number of cells
        unsigned nel = solid_mesh_pt()->nelement();

        // set up vtu Points
        struct Point
        {
            Vector<double> coordinate;
            struct Data
            {
                std::string name;
                std::vector<double> value;
            };
            std::vector<Data> data;
        };
        std::vector<Point> points;
        points.reserve(nel * sList.size());

        // set up vtu Cells
        struct Cell
        {
            std::vector<unsigned long> connectivity;
            unsigned long offset;
            unsigned type;
        };
        std::vector<Cell> cells;
        points.reserve(nel);

        //for all elements
        for (unsigned e = 0; e < nel; e++)
        {
            // Get pointer to element
            auto el_pt = dynamic_cast<ELEMENT*>(
                solid_mesh_pt()->element_pt(e));

            std::vector<unsigned long> connectivity;
            unsigned ix = 0; // node index
            for (const auto& s : sList)
            {
                // pointer to node info
                auto n = nList[ix];
                auto n_pt = dynamic_cast<SolidNode*>(el_pt->node_pt(n));
                // Get Eulerian coordinates
                Vector<double> x(3);
                el_pt->interpolated_x(s, x);
                // get velocity
                Vector<double> dxdt(3);
                el_pt->interpolated_dxdt(s, 1, dxdt);
                // get displacement
                Vector<double> xi(3);
                el_pt->interpolated_xi(s, xi);
                // getBodyForce
                Vector<double> body_force(3);
                auto bodyForceFct = el_pt->body_force_fct_pt();
                if (bodyForceFct) bodyForceFct(time(), xi, body_force);
                // get stress/pressure
                DenseMatrix<double> sigma(3,3);
                el_pt->get_stress(s, sigma);
                double pressure = (sigma(0,0)+sigma(1,1)+sigma(2,2))/3;
                // get strain
                DenseMatrix<double> strain(3,3);
                el_pt->get_strain(s, strain);
                // get velocity (fails)
                std::vector<double> dudt
                    {el_pt->dnodal_position_dt(n, 0),
                     el_pt->dnodal_position_dt(n, 1),
                     el_pt->dnodal_position_dt(n, 2)};
                // get boundary (works)
                std::set<unsigned>* boundaries_pt;
                n_pt->get_boundaries_pt(boundaries_pt);
                double b = boundaries_pt ? *boundaries_pt->begin() : -1;
                std::vector<double> pin {(double) n_pt->position_is_pinned(0),
                                         (double) n_pt->position_is_pinned(1),
                                         (double) n_pt->position_is_pinned(2)};
                points.push_back(
                    {x, {
                             {"Velocity", {dxdt[0], dxdt[1], dxdt[2]}},
                             {"Displacement", {x[0] - xi[0], x[1] - xi[1], x[2] - xi[2]}},
                             {"BodyForce", {body_force[0], body_force[1], body_force[2]}},
                             {"Pin", pin},
                             {"Velocity2", dudt},
//                             {"Pressure", {pressure}},
//                             {"Stress", {sigma(0,0), sigma(0,1), sigma(0,2),
//                                           sigma(1,0), sigma(1,1), sigma(1,2),
//                                           sigma(2,0), sigma(2,1), sigma(2,2)}},
//                             {"Strain", {strain(0,0), strain(0,1), strain(0,2),
//                                           strain(1,0), strain(1,1), strain(1,2),
//                                           strain(2,0), strain(2,1), strain(2,2)}},
                             //\todo undo this comment, but it causes a problem in the output
                             {"Boundary", {b}}
                         }});

                // add to connectivity
                connectivity.push_back(points.size() - 1);
                ix++;
            }
            cells.push_back({connectivity, points.size(), vtkFormat});
        }

        //open vtk file
        static unsigned count = 0;
#ifdef OOMPH_HAS_MPI
        std::string vtkFileName = name_ + "FEM_" + std::to_string(OOMPH_MPI_PROCESSOR_ID) + "_" + std::to_string(count++) + ".vtu";
#else
        std::string vtkFileName = name_ + "FEM_" + std::to_string(count++) + ".vtu";
#endif

        std::ofstream vtk(vtkFileName);

        //write vtk file
        vtk << "<?xml version=\"1.0\"?>\n"
               "<!-- time 10.548-->\n"
               "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
               "<UnstructuredGrid>\n"
               "<Piece NumberOfPoints=\""
            << points.size()
            << "\" NumberOfCells=\""
            << cells.size()
            << "\">\n"
               "\n"
               "<Points>\n"
               "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (auto point : points)
        {
            vtk << point.coordinate[0] << " " << point.coordinate[1] << " " << point.coordinate[2] << "\n";
        }
        vtk << "</DataArray>\n"
               "</Points>\n\n";
        if (not points.empty())
        {
            vtk << "<PointData  Vectors=\"vector\">\n";
            for (int i = 0; i < points.front().data.size(); ++i)
            {
                auto data = points.front().data[i];
                vtk << R"(<DataArray type="Float32" Name=")"
                    << points.front().data[i].name
                    << R"(" NumberOfComponents=")"
                    << points.front().data[i].value.size()
                    << R"(" format="ascii">)" << "\n";
                for (const Point& point : points)
                {
                    for (const auto& value : point.data[i].value)
                    {
                        vtk << value << " ";
                    }
                }
                vtk << "\n"
                       "</DataArray>\n";
            }
            vtk << "</PointData>\n"
                   "\n";
        }
        vtk << "<Cells>\n"
               "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (const Cell& cell : cells)
        {
            for (auto point : cell.connectivity)
            {
                vtk << point << " ";
            }
        }
        vtk << "</DataArray>\n"
               "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (const Cell& cell : cells)
        {
            vtk << cell.offset << " ";
        }
        vtk << "</DataArray>\n"
               "<DataArray type=\"UInt8\"  Name=\"types\" format=\"ascii\">\n";
        for (const Cell& cell : cells)
        {
            vtk << cell.type << " ";
        }
        vtk << "\n"
               "</DataArray>\n"
               "</Cells>\n"
               "\n"
               "</Piece>\n"
               "</UnstructuredGrid>\n"
               "</VTKFile>\n";
        vtk.close();
        logger(INFO, "Written %", vtkFileName);
    }
};
#endif //SOLID_PROBLEM_H
