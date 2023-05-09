//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef VOLUME_COUPLING_H
#define VOLUME_COUPLING_H

#include "BaseCoupling.h"

class VolumeCoupling : public BaseCoupling
{
public:
    typedef VolumeCoupledElement<OomphProblem::ELEMENT> VELEMENT;

    /**
     * information needed to be stored to link surface-coupled elements to DPM
     */
    struct DPMVCoupledElement
    {
        // the pointer to a bulk element
        VELEMENT* bulk_elem_pt;
        // the vector that contains pointers to the coupled particles
        Vector<BaseParticle*> listOfCoupledParticles;
        Vector<BaseParticle*> listOfCoupledParticlesExt;
        // the list of local coordinates of particle centers
        Vector<Vector<double>> listOfParticleCentersLoc;
        // the list of shape functions evaluated at particle centers
        Vector<Vector<double>> listOfShapesAtParticleCenters;
    };

    // default constructor
    VolumeCoupling() = default;

    // constructor passing in the minimal coupling weight and flag for using the implicit or explicit scheme
    VolumeCoupling(const double& wMin, const bool& flag)
    {
        setExplicit(flag);
        weightMin_ = wMin;
    }

    /**
 *  Override computeExternalForces to add the vcoupling forces
 */
    void computeExternalForces(BaseParticle* CI) override
    {
        //Checks that the current particle is not "fixed"
        //and hence infinitely massive!
        if (!CI->isFixed())
        {
            // Applying the force due to bodyForce_fct (F = m.g)
            CI->addForce(getGravity() * CI->getMass());
            // add particle coupling force 1/w_i*f_i^c to the total force
            CI->addForce(CI->getCouplingForce());
        }
    }

    /**
     * Override actionsAfterTimeStep to reset vcoupling force
     * \todo actions after time step should be user-defined, so not be used here
     */
    void actionsAfterTimeStep() override
    {
        // Reset coupling weights and forces to default values
        for (BaseParticle* p : particleHandler)
        {
            if (!p->isFixed())
            {
                // reset the coupling force to zero for the next time step
                p->resetCouplingForce();
                p->resetInvCouplingWeight();
                p->resetCouplingMass();
            }
        }

        // check how much equilibrium is broken
        const auto p = particleHandler.getLastObject();
        Vec3D f = p->getMass() * getGravity();
        Vec3D a = p->getForce();
        if (a.getLengthSquared() > 1.0e-5 * f.getLengthSquared())
        {
            logger(INFO, "time %, body force: %, total force: %, position %",
                   getTime(), f, a, p->getPosition());
        }

        // check if each particle-wall interaction is unique
        for (auto p : particleHandler)
        {
            unsigned n = 0;
            for (auto i : p->getInteractions())
            {
                if (i->getI()->getName() == "TriangleWall") n++;
            }
            if (n > 1)
            {
                logger(INFO, "particle % is interacting with % TriangleWalls", p->getId(), n);
                for (auto i : p->getInteractions())
                {
                    if (i->getI()->getName() == "TriangleWall")
                    {
                        logger(INFO, "particle % is interacting with wall %", p->getId(), i->getI()->getId());
                    }
                }
            }
        }
    }

    // solve volume coupled OomphMercuryProblem using CG for micro-macro mapping
    void solveVolumeCoupling(const double& width)
    {
        // set the coarse-grainning width w.r.t. length scale
        setCGWidth(width * particleHandler.getLargestParticle()->getRadius());
        // solve volume coupled OomphMercuryProblem
        solveVolumeCoupling();
    }

    // solve volume coupled OomphMercuryProblem
    void solveVolumeCoupling()
    {
        initialiseSolve();

        // get the time step of mercuryProb
        double mercury_dt = getTimeStep();
        /*
                // %TODO use the DPM time scale to non-dimensionalise the OomphProblem<SELEMENT,TIMESTEPPER>?
                Global_Physical_Variables::timeScale = mercury_dt;
        */

        /*      // get number of mercury-steps per Oomph-step (critical time step = intrinsic_time or from_eigen)
                double oomph_dt = 0.2*Global_Physical_Variables::intrinsic_time;
                double oomph_dt = 0.2*OomphProblem::get_critical_timestep_from_eigen();
                unsigned nstep = unsigned(floor(oomph_dt/mercury_dt));
        */
        unsigned nstep = 1;
        double oomph_dt = nstep * mercury_dt;

        // setup initial conditions
        OomphProblem::set_initial_conditions(oomph_dt);
        //        MercuryBeam::set_initial_conditions();

        // This is the main loop over advancing time
        unsigned nDone = 0;
        while (getTime() < getTimeMax())
        {
            // solve the coupled problem for one time step
            computeOneTimeStepForVCoupling(nstep);
            // write outputs of the oomphProb
            if (getVtkWriter()->getFileCounter() > nDone)
            {
                OomphProblem::doc_solution();
                nDone = getVtkWriter()->getFileCounter();
            }
        }
        // close output files of mercuryProb
        finaliseSolve();
    }

    // solve volume coupled OomphMercuryProblem for one time step
    void computeOneTimeStepForVCoupling(const unsigned& nstep)
    {
        listOfDPMVCoupledElements_.clear();
        checkParticlesInFiniteElems();
        //        setNodalWeightOnCubicMesh();
        setNodalWeightByDistance();
        computeWeightOnParticles();
        if (explicit_)
        {
            computeCouplingForce();
            solve();
            solveMercury(nstep);
        }
        else
        {
            solveFirstHalfTimeStep();
            computeCouplingOnFEM();
            solveOomph();
            computeCouplingOnDPM();
            solveSecondHalfTimeStep();
            solveMercury(nstep - 1);
        }
    }

    /**
     * Used in OomphMercuryCoupling::computeOneTimeStepForVCoupling to do the actions before time step
     */
    void solveFirstHalfTimeStep()
    {
        logger(DEBUG, "starting solveHalfTimeStep()");

        logger(DEBUG, "about to call writeOutputFiles()");
        writeOutputFiles(); //everything is written at the beginning of the time step!

        logger(DEBUG, "about to call hGridActionsBeforeIntegration()");
        hGridActionsBeforeIntegration();

        //Computes the half-time step velocity and full time step position and updates the particles accordingly
        logger(DEBUG, "about to call integrateBeforeForceComputation()");
        integrateBeforeForceComputation();
        //New positions require the MPI and parallel periodic boundaries to do things
        logger(DEBUG, "about to call performGhostParticleUpdate()");
        performGhostParticleUpdate();

        /// \todo MX: this is not true anymore. all boundaries are handled here.
        /// particles have received a position update, so here the deletion boundary deletes particles
        ///\TODO add particles need a periodic check

        logger(DEBUG, "about to call checkInteractionWithBoundaries()");
        checkInteractionWithBoundaries(); // INSERTION boundaries handled

        logger(DEBUG, "about to call hGridActionsAfterIntegration()");
        hGridActionsAfterIntegration();

        // Compute forces
        ///\bug{In chute particles are added in actions_before_time_set(), however they are not written to the xballs data yet, but can have a collision and be written to the fstat data}
        // INSERTION/DELETION boundary flag change
        for (BaseBoundary* b : boundaryHandler)
        {
            b->checkBoundaryBeforeTimeStep(this);
        }

        logger(DEBUG, "about to call actionsBeforeTimeStep()");
        actionsBeforeTimeStep();

        logger(DEBUG, "about to call checkAndDuplicatePeriodicParticles()");
        checkAndDuplicatePeriodicParticles();

        logger(DEBUG, "about to call hGridActionsBeforeTimeStep()");
        hGridActionsBeforeTimeStep();
    }

    /**
     * Used in OomphMercuryCoupling::computeOneTimeStepForVCoupling to do the actions after time step
     */
    void solveSecondHalfTimeStep()
    {
        //Creates and updates interactions and computes forces based on these
        logger(DEBUG, "about to call computeAllForces()");
        computeAllForces();

        logger(DEBUG, "about to call removeDuplicatePeriodicParticles()");
        removeDuplicatePeriodicParticles();

        //Computes new velocities and updates the particles accordingly
        logger(DEBUG, "about to call integrateAfterForceComputation()");
        integrateAfterForceComputation();

        logger(DEBUG, "about to call actionsAfterTimeStep()");
        actionsAfterTimeStep();

        //erase interactions that have not been used during the last time step
        logger(DEBUG, "about to call interactionHandler.eraseOldInteractions(getNumberOfTimeSteps())");
        interactionHandler.eraseOldInteractions(getNumberOfTimeSteps());
        logger(DEBUG, "about to call interactionHandler.actionsAfterTimeStep()");
        interactionHandler.actionsAfterTimeStep();
        particleHandler.actionsAfterTimeStep();

        continueInTime();

        logger(DEBUG, "finished computeOneTimeStep()");
    }

    void checkParticlesInFiniteElems()
    {
        // how many bulk elements in total
        unsigned n_element = OomphProblem::solid_mesh_pt()->nelement();

        // loop over the bulk elements
        for (unsigned e = 0; e < n_element; e++)
        {
            // get pointer to the bulk element
            auto bulk_elem_pt = dynamic_cast<VELEMENT*>(OomphProblem::solid_mesh_pt()->element_pt(e));

            // get number of nodes in the bulk element
            const unsigned nnode = bulk_elem_pt->nnode();

            // Set up memory for the shape/test functions
            Shape psi(nnode);

            // three arrays that contain the x, y and z coordinates of the bulk element
            Vector<double> listOfCoordX;
            Vector<double> listOfCoordY;
            Vector<double> listOfCoordZ;

            // get the x, y and z coordinates of the bulk element
            for (unsigned n = 0; n < nnode; n++)
            {
                listOfCoordX.push_back(bulk_elem_pt->node_pt(n)->x(0) * Global_Physical_Variables::lenScale);
                listOfCoordY.push_back(bulk_elem_pt->node_pt(n)->x(1) * Global_Physical_Variables::lenScale);
                listOfCoordZ.push_back(bulk_elem_pt->node_pt(n)->x(2) * Global_Physical_Variables::lenScale);
            }

            // get the bounding box of the bulk element
            Vec3D min;
            min.X = *min_element(listOfCoordX.begin(), listOfCoordX.end());
            min.Y = *min_element(listOfCoordY.begin(), listOfCoordY.end());
            min.Z = *min_element(listOfCoordZ.begin(), listOfCoordZ.end());
            Vec3D max;
            max.X = *max_element(listOfCoordX.begin(), listOfCoordX.end());
            max.Y = *max_element(listOfCoordY.begin(), listOfCoordY.end());
            max.Z = *max_element(listOfCoordZ.begin(), listOfCoordZ.end());

            // extend the bounding box if construct mapping with coarse graining
            min -= Vec3D(1.0, 1.0, 1.0) * ( getCGWidth() - 2 * particleHandler.getLargestParticle()->getRadius() );
            max += Vec3D(1.0, 1.0, 1.0) * ( getCGWidth() - 2 * particleHandler.getLargestParticle()->getRadius() );

            // get particles if there are in the bounding box
            Vector<BaseParticle*> pList;
            getParticlesInCell(min, max, pList);
            if (pList.size() != 0)
            {
                DPMVCoupledElement DPM_coupled_elem;
                // save the extended particle list if construct mapping with coarse graining functions
                if (useCGMapping())
                { DPM_coupled_elem.listOfCoupledParticlesExt = pList; }
                // check if a particle resides within the bulk element
                for (auto it = pList.begin(); it != pList.end(); it++)
                {
                    Vector<double> s(3, 0.0), x(3, 0.0);
                    GeomObject* geom_obj_pt = 0;
                    // get global coordinates of particle center (DPM)
                    Vec3D xp = ( *it )->getPosition();
                    xp /= Global_Physical_Variables::lenScale;
                    x[0] = xp.getX();
                    x[1] = xp.getY();
                    x[2] = xp.getZ();
                    // get local coordinates of particle center (FEM)
                    bulk_elem_pt->locate_zeta(x, geom_obj_pt, s);
                    // If true, the iterator is decremented after it is passed to erase() but before erase() is executed
                    if (geom_obj_pt == nullptr)
                    { pList.erase(it--); }
                    else
                    {
                        // save local coordinates of particle center
                        DPM_coupled_elem.listOfParticleCentersLoc.push_back(s);
                        // evaluate shape functions at particle center
                        bulk_elem_pt->shape(s, psi);
                        Vector<double> shapes;
                        for (unsigned l = 0; l < nnode; l++) shapes.push_back(psi(l));
                        // save shape functions at particle center
                        DPM_coupled_elem.listOfShapesAtParticleCenters.push_back(shapes);
                    }
                }
                // if there are particles in the bulk element (check with locate_zeta)
                if (pList.size() != 0)
                {
                    // save pointer to the bulk element and the pointer list to particles
                    DPM_coupled_elem.bulk_elem_pt = bulk_elem_pt;
                    DPM_coupled_elem.listOfCoupledParticles = pList;
                    listOfDPMVCoupledElements_.push_back(DPM_coupled_elem);
                }
            }
        }
    }

    /**
     * VolumeCoupling: setNodalWeightOnCubicMesh
     */
    void setNodalWeightOnCubicMesh()
    {
        // get a full list of all nodes that are coupled
        Vector<CoupledSolidNode*> listOfSolidNodes;
        for (DPMVCoupledElement coupled_elem : listOfDPMVCoupledElements_)
        {
            // get pointers to all nodes of the bulk element
            const unsigned nnode = coupled_elem.bulk_elem_pt->nnode();
            for (unsigned n = 0; n < nnode; n++)
            {
                CoupledSolidNode* node_pt = dynamic_cast<CoupledSolidNode*>(coupled_elem.bulk_elem_pt->node_pt(n));
                // if it is not a boundary node, add it to the list
                if (!node_pt->is_on_boundary())
                {
                    listOfSolidNodes.push_back(node_pt);
                }
                    // set nodal weights to one if it is a boundary node
                else
                { node_pt->set_coupling_weight(0.5); }
            }
        }
        // count number of occurrences of a node and set nodal weights
        for (CoupledSolidNode* node_pt : listOfSolidNodes)
        {
            unsigned n = count(listOfSolidNodes.begin(), listOfSolidNodes.end(), node_pt);
            if (n == pow(2, getSystemDimensions()))
            { node_pt->set_coupling_weight(0.5); }
            else
            { node_pt->set_coupling_weight(0.5); }
            //            cout << "The solid node " << node_pt << " is shared by " << n << " coupled elements. ";
            //            cout << "Set weight to " << node_pt->get_coupling_weight() << endl;
        }
    }

    /**
     * used in computeOneTimeStepForVCoupling (VolumeCoupling)
     */
    void setNodalWeightByDistance()
    {
        unsigned direction = 1;
        // store the position[direction] of all coupled nodes in a vector
        Vector<double> listOfPos;
        for (DPMVCoupledElement coupled_elem : listOfDPMVCoupledElements_)
        {
            // get pointers to all nodes of the bulk element
            const unsigned nnode = coupled_elem.bulk_elem_pt->nnode();
            for (unsigned n = 0; n < nnode; n++)
            {
                CoupledSolidNode* node_pt = dynamic_cast<CoupledSolidNode*>(coupled_elem.bulk_elem_pt->node_pt(n));
                listOfPos.push_back(node_pt->x(direction));
            }
        }
        // get the minimum and maximum coordinates[direction]
        auto posMin = *min_element(listOfPos.begin(), listOfPos.end());
        auto posMax = *max_element(listOfPos.begin(), listOfPos.end());
        /*!
         * assign coupling weight to coupled solid nodes
         * implicit volume coupling allows the minimal coupling weight to be very small
         * e.g. double weightMin = 0.1 for the explicit volume coupling scheme
         */
        double weightMax = 1.0 - weightMin_;
        for (DPMVCoupledElement coupled_elem : listOfDPMVCoupledElements_)
        {
            // get pointers to all nodes of the bulk element
            const unsigned nnode = coupled_elem.bulk_elem_pt->nnode();
            for (unsigned n = 0; n < nnode; n++)
            {
                CoupledSolidNode* node_pt = dynamic_cast<CoupledSolidNode*>(coupled_elem.bulk_elem_pt->node_pt(n));
                // assign nodal coupling weight using a linear function of distance
                double weight = weightMin_
                    + ( weightMax - weightMin_ ) / ( posMax - posMin ) * ( node_pt->x(direction) - posMin );
                //                // assign nodal coupling weight using a cosine function of distance
                //                double amp = 0.5*(weightMax - weightMin_);
                //                double omega = M_PI/(posMax-posMin);
                //                double weight = weightMin_ + amp - amp*cos(omega*(node_pt->x(direction)-posMin));
                node_pt->set_coupling_weight(weight);
                //                cout << "Set weight to " << node_pt->get_coupling_weight();
                //                cout << " at a distance of " <<  node_pt->x(direction)-posMin <<endl;
            }
        }
    }

    void getProjectionAndProjected(const DPMVCoupledElement& coupled_elem,
                                   Vector<Vector<double>>& couplingMatrix, Vector<Vector<double>>& nodalDisplDPM)
    {
        // get coupled particles in the bulk element
        Vector<BaseParticle*> particles;
        if (!useCGMapping())
        { particles = coupled_elem.listOfCoupledParticles; }
        else
        { particles = coupled_elem.listOfCoupledParticlesExt; }
        // get number of coupled particles in the bulk element
        const unsigned nParticles = particles.size();
        // get nodal coordinates of the bulk element
        const unsigned nnode = coupled_elem.bulk_elem_pt->nnode();
        // get dimension of the problem
        const unsigned dim = coupled_elem.bulk_elem_pt->dim();
        // prepare vector of particle displacement
        Vector<Vector<double>> particleDispl(nParticles, Vector<double>(dim, 0));
        // fill in the projection matrix of the coupled bulk element
        getProjectionMatrix(coupled_elem, couplingMatrix);
        // get particle displacements
        for (unsigned m = 0; m < nParticles; m++)
        {
            Vec3D displ = particles[m]->getDisplacement() / Global_Physical_Variables::lenScale;
            particleDispl[m][0] = displ.getX();
            particleDispl[m][1] = displ.getY();
            particleDispl[m][2] = displ.getZ();
        }
        // get projected nodal displacement field from particle displacement at the center
        for (unsigned i = 0; i < nnode; i++)
        {
            for (unsigned k = 0; k < dim; k++)
            {
                for (unsigned j = 0; j < nParticles; j++)
                {
                    nodalDisplDPM[i][k] += couplingMatrix[i][j] * particleDispl[j][k];
                }
            }
        }
    }

    void getProjectionMatrix(const DPMVCoupledElement& coupled_elem, Vector<Vector<double>>& couplingMatrix)
    {
        // get number of coupled particles in the bulk element
        const unsigned nParticles = couplingMatrix[0].size();
        // get nodal coordinates of the bulk element
        const unsigned nnode = couplingMatrix.size();
        const unsigned dim = coupled_elem.bulk_elem_pt->dim();
        // if construct mapping with FEM basis functions
        if (!useCGMapping())
        {
            // loop over shape functions at the nodes
            for (unsigned l = 0; l < nnode; l++)
            {
                double sum = 0;
                // loop over shape functions at the particles
                for (unsigned m = 0; m < nParticles; m++)
                {
                    // get the l-th shape function evaluated at the center of particle m
                    double shape = coupled_elem.listOfShapesAtParticleCenters[m][l];
                    // shape function weighted by particle volume N_{l,m}*V_m
                    shape *= coupled_elem.listOfCoupledParticles[m]->getVolume();
                    couplingMatrix[l][m] = shape;
                    sum += shape;
                }
                // normalize each row of the projection matrix to sum to one
                for (unsigned m = 0; m < nParticles; m++)
                {
                    // projection rule defined as N_{l,m}*V_m / sum_m N_{l,m}*V_m
                    if (sum != 0)
                    { couplingMatrix[l][m] /= sum; }
                }
            }
        }
        else
        {
            // set up memory for the shape functions
            Shape psi(nnode);
            // set up memory for the local and global coordinates
            Vector<double> s(dim), x(dim);
            // get the element pointer
            auto el_pt = coupled_elem.bulk_elem_pt;
            // get the number of integration points
            const unsigned n_intpt = el_pt->integral_pt()->nweight();
            // loop over the integration points
            for (unsigned ipt = 0; ipt < n_intpt; ipt++)
            {
                // get the integral weight
                double w = el_pt->integral_pt()->weight(ipt);
                // find the shape functions at the integration points r_i
                el_pt->shape_at_knot(ipt, psi);
                // get the value of the local coordinates at the integration point
                for (unsigned i = 0; i < dim; i++)
                { s[i] = el_pt->integral_pt()->knot(ipt, i); }
                // get the value of the global coordinates at the integration point
                el_pt->interpolated_x(s, x);
                // loop over the nodes
                for (unsigned l = 0; l < nnode; l++)
                {
                    // loop over shape functions at the particles
                    for (unsigned m = 0; m < nParticles; m++)
                    {
                        // set CG coordinates
                        CGCoordinates::XYZ coord;
                        coord.setXYZ(Vec3D(x[0], x[1], x[2]));
                        // evaluate the value of CG function around particle m at CGcoords \phi(\vec r_i-r_m)
                        double phi = getCGFunction().evaluateCGFunction(
                            coupled_elem.listOfCoupledParticlesExt[m]->getPosition() /
                                Global_Physical_Variables::lenScale, coord);
                        // CG mapping defined as \tilde{N_{l,m}}_ipt = w_ipt * \phi(\vec r_i-r_m) * N_l(r_i)
                        double shape = w * phi * psi(l);
                        // CG mapping weighted by particle volume \tilde{N_{l,m}}*V_m
                        shape *= coupled_elem.listOfCoupledParticlesExt[m]->getVolume();
                        // sum over the integration points
                        couplingMatrix[l][m] += shape;
                    }
                }
            }
            // normalize each row of the projection matrix to sum to one
            for (unsigned l = 0; l < nnode; l++)
            {
                double sum = 0;
                for (unsigned m = 0; m < nParticles; m++)
                { sum += couplingMatrix[l][m]; }
                for (unsigned m = 0; m < nParticles; m++)
                {
                    // projection rule defined as \tilde{N_{l,m}}*V_m / sum_m \tilde{N_{l,m}}*V_m
                    if (sum != 0)
                    { couplingMatrix[l][m] /= sum; }
                }
            }
        }
    }

    void computeCouplingForce()
    {
        // first loop over the coupled bulk elements
        for (DPMVCoupledElement coupled_elem : listOfDPMVCoupledElements_)
        {
            // get pointer to the bulk element
            VELEMENT* elem_pt = coupled_elem.bulk_elem_pt;
            // get particles and the number of coupled particles in the bulk element
            Vector<BaseParticle*> particles;
            if (!useCGMapping())
            { particles = coupled_elem.listOfCoupledParticles; }
            else
            { particles = coupled_elem.listOfCoupledParticlesExt; }
            const unsigned nParticles = particles.size();
            // get number of nodes in the bulk element
            const unsigned nnode = coupled_elem.bulk_elem_pt->nnode();
            // prepare projection matrix
            Vector<Vector<double>> couplingMatrix(nnode, Vector<double>(nParticles, 0.0));
            // prepare vector of nodal displacements projected from particle centers
            const unsigned dim = elem_pt->dim();
            Vector<Vector<double>> nodalDisplDPM(nnode, Vector<double>(dim, 0.0));
            // compute projected nodal displacements
            getProjectionAndProjected(coupled_elem, couplingMatrix, nodalDisplDPM);
            // compute nodal coupling forces by penalizing the difference in displacement
            Vector<Vector<double>> nodalCouplingForces(nnode, Vector<double>(dim, 0.0));
            computeNodalCouplingForces(elem_pt, nnode, dim, nodalDisplDPM, nodalCouplingForces);
            // get coupling forces projected from coupled nodes to particles (note couplingMatrix^{-1} = couplingMatrix^T)
            Vector<Vector<double>> particleCouplingForces(nParticles, Vector<double>(dim, 0.0));
            for (unsigned i = 0; i < nParticles; i++)
            {
                for (unsigned k = 0; k < dim; k++)
                {
                    for (unsigned j = 0; j < nnode; j++)
                    {
                        // note that particle and nodal coupling forces are in opposite directions
                        // the negative sign on the FEM side is taken care of in by element
                        particleCouplingForces[i][k] += couplingMatrix[j][i] * nodalCouplingForces[j][k];
                    }
                }
                // add particle coupling forces to be used by actionsAfterTimeStep(...)
                auto p = particles[i];
                Vec3D cForce = Vec3D(particleCouplingForces[i][0],
                                     particleCouplingForces[i][1],
                                     particleCouplingForces[i][2]);
                // if construct mapping with coarse graining, add elements' contribution to particle coupling forces
                if (!useCGMapping())
                {
                    p->setCouplingForce(cForce * Global_Physical_Variables::forceScale());
                }
                    // otherwise, set particle coupling force since the operation is locally within each element
                else
                {
                    p->addCouplingForce(cForce * Global_Physical_Variables::forceScale());
                }
            }
            // assign nodal coupling force to the element to be used by element::fill_in_contribution_to_residuals(...)
            elem_pt->set_nodal_coupling_residual(nodalCouplingForces);
        }
    }

    void computeCouplingOnFEM()
    {
        // first loop over the coupled bulk elements
        for (DPMVCoupledElement coupled_elem : listOfDPMVCoupledElements_)
        {
            // get pointer to the bulk element
            VELEMENT* elem_pt = coupled_elem.bulk_elem_pt;

            // get particles and the number of coupled particles in the bulk element
            Vector<BaseParticle*> particles;
            if (!useCGMapping())
            { particles = coupled_elem.listOfCoupledParticles; }
            else
            { particles = coupled_elem.listOfCoupledParticlesExt; }
            const unsigned nParticles = particles.size();

            // get number of nodes in the bulk element
            const unsigned nnode = coupled_elem.bulk_elem_pt->nnode();
            const unsigned dim = elem_pt->dim();

            // prepare projection matrix \Pi_{i,\alpha}
            Vector<Vector<double>> couplingMatrix(nnode, Vector<double>(nParticles, 0.0));
            // fill in the projection matrix of the coupled bulk element
            getProjectionMatrix(coupled_elem, couplingMatrix);

            // prepare the coupled stiffness matrix of the bulk element, penalty*V_{ij}
            Vector<Vector<double>> bulkStiffness(nnode, Vector<double>(nnode, 0.0));
            // fill in the volume coupled stiffness matrix
            getBulkStiffness(elem_pt, bulkStiffness);
            // assign penalty*\Pi_{\alpha,i}*V_{i,j} to the VCoupled element
            Vector<Vector<double>> couplingStiffness(nParticles, Vector<double>(nnode, 0.0));
            for (unsigned m = 0; m < nParticles; m++)
            {
                for (unsigned j = 0; j < nnode; j++)
                {
                    for (unsigned i = 0; i < nnode; i++)
                    {
                        couplingStiffness[m][j] += couplingMatrix[i][m] * bulkStiffness[i][j];
                    }
                }
            }
            elem_pt->setCouplingStiffness(couplingStiffness);

            // prepare the coupling residual projected from the particles to the node
            Vector<Vector<double>> nodalResidualDPM(nnode, Vector<double>(dim, 0.0));

            // get displacement u(t+dt) = v(t+0.5dt)*dt after integrateBeforeForceComputation(...)
            for (unsigned m = 0; m < nParticles; m++)
            {
                Vec3D cForce;
                for (unsigned l = 0; l < nnode; l++)
                {
                    for (unsigned ll = 0; ll < nnode; ll++)
                    {
                        for (unsigned mm = 0; mm < nParticles; mm++)
                        {
                            // change finite element coupling stiffness into the DEM scale
                            cForce += couplingMatrix[l][m] * bulkStiffness[l][ll] *
                                Global_Physical_Variables::stiffScale() * couplingMatrix[ll][mm] *
                                particles[mm]->getDisplacement();
                        }
                    }
                }
                // add the first particle coupling force: -penalty*\PI_{\alpha i}*V_{ij}*\PI_{j \beta}*u_\beta
                particles[m]->addCouplingForce(-cForce);
            }
            // compute the residual projected from the particles to the nodes
            for (unsigned i = 0; i < nnode; i++)
            {
                for (unsigned j = 0; j < nnode; j++)
                {
                    for (unsigned m = 0; m < nParticles; m++)
                    {
                        Vec3D displ = particles[m]->getDisplacement() / Global_Physical_Variables::lenScale;
                        double particleStiffness = bulkStiffness[i][j] * couplingMatrix[j][m];
                        nodalResidualDPM[i][0] -= particleStiffness * displ.getX();
                        nodalResidualDPM[i][1] -= particleStiffness * displ.getY();
                        nodalResidualDPM[i][2] -= particleStiffness * displ.getZ();
                    }
                }
            }
            // assign residual to the coupled bulk element to be used by element::fill_in_contribution_to_residuals(...)
            elem_pt->set_nodal_coupling_residual(nodalResidualDPM);

            // assign Jacobian to the coupled bulk element to be used by element::fill_in_contribution_to_residuals(...)
            elem_pt->set_nodal_coupling_jacobian(bulkStiffness);
        }
    }

    void computeCouplingOnDPM()
    {
        // first loop over the coupled bulk elements
        for (DPMVCoupledElement coupled_elem : listOfDPMVCoupledElements_)
        {
            // get pointer to the bulk element
            VELEMENT* elem_pt = coupled_elem.bulk_elem_pt;

            // get particles and the number of coupled particles in the bulk element
            Vector<BaseParticle*> particles;
            if (!useCGMapping())
            { particles = coupled_elem.listOfCoupledParticles; }
            else
            { particles = coupled_elem.listOfCoupledParticlesExt; }
            const unsigned nParticles = particles.size();

            // get number of nodes in the bulk element
            const unsigned nnode = coupled_elem.bulk_elem_pt->nnode();
            const unsigned dim = elem_pt->dim();

            // loop over the particles
            for (unsigned m = 0; m < nParticles; m++)
            {
                // get "residual" from the bulk element
                Vec3D cForce;
                for (unsigned l = 0; l < nnode; l++)
                {
                    double coupling = elem_pt->getCouplingStiffness(m, l);
                    // get the nodal displacement
                    Vec3D force = coupling * Vec3D(elem_pt->nodal_position(0, l, 0) - elem_pt->nodal_position(1, l, 0),
                                                   elem_pt->nodal_position(0, l, 1) - elem_pt->nodal_position(1, l, 1),
                                                   elem_pt->nodal_position(0, l, 2) - elem_pt->nodal_position(1, l, 2));
                    // add contributions to the coupling force
                    cForce += force;
                }
                // add the second particle coupling force: penalty*\PI_{\alpha i}*V_{ij}*u_j
                particles[m]->addCouplingForce(cForce * Global_Physical_Variables::forceScale());
            }
        }
    }

    void computeNodalCouplingForces(VELEMENT*& elem_pt, const unsigned& nnode, const unsigned& dim,
                                    const Vector<Vector<double>>& nodalDisplDPM,
                                    Vector<Vector<double>>& nodalCouplingForces)
    {
        // get the difference between nodal displacement from FEM and projected ones from DPM
        Vector<Vector<double>> nodalDisplDiff(nnode, Vector<double>(dim, 0));
        for (unsigned l = 0; l < nnode; l++)
        {
            for (unsigned i = 0; i < dim; i++)
            {
                double displFEM = elem_pt->nodal_position(0, l, i) - elem_pt->nodal_position(1, l, i);
                double displDPM = nodalDisplDPM[l][i];
                nodalDisplDiff[l][i] = displDPM - displFEM;
            }
        }
        // Set up memory for the shape functions
        Shape psi(nnode);
        DShape dpsidxi(nnode, dim);
        const unsigned n_intpt = elem_pt->integral_pt()->nweight();
        // Loop over the integration points
        for (unsigned ipt = 0; ipt < n_intpt; ipt++)
        {
            // Get the integral weight
            double w = elem_pt->integral_pt()->weight(ipt);
            // Call the derivatives of the shape functions (and get Jacobian)
            double J = elem_pt->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);
            // Loop over the test functions, nodes of the element
            for (unsigned l = 0; l < nnode; l++)
            {
                // Loop over the nodes of the element again
                for (unsigned ll = 0; ll < nnode; ll++)
                {
                    double vol = Global_Physical_Variables::penalty * psi(ll) * psi(l) * w * J;
                    /*
                                        for (unsigned i=0; i<dim; i++)
                                        {
                                            vol += Global_Physical_Variables::penalty*2*Global_Physical_Variables::L/3.0 * dpsidxi(ll,i)*dpsidxi(l,i)*w*J;
                                        }
                    */
                    for (unsigned i = 0; i < dim; i++)
                    {
                        nodalCouplingForces[l][i] += -vol * nodalDisplDiff[l][i];
                    }
                }
            }
        }
    }


    void getBulkStiffness(VELEMENT*& elem_pt, Vector<Vector<double>>& bulkStiffness)
    {
        // get number of nodes in the bulk element
        const unsigned nnode = elem_pt->nnode();
        // get dimension of the problem
        const unsigned dim = elem_pt->dim();
        // Set up memory for the shape functions
        Shape psi(nnode);
        DShape dpsidxi(nnode, dim);
        const unsigned n_intpt = elem_pt->integral_pt()->nweight();
        // Loop over the integration points
        for (unsigned ipt = 0; ipt < n_intpt; ipt++)
        {
            // Get the integral weight
            double w = elem_pt->integral_pt()->weight(ipt);
            // Call the derivatives of the shape functions (and get Jacobian)
            double J = elem_pt->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);
            // Loop over the test functions, nodes of the element
            for (unsigned l = 0; l < nnode; l++)
            {
                // Loop over the nodes of the element again
                for (unsigned ll = 0; ll < nnode; ll++)
                {
                    bulkStiffness[l][ll] += Global_Physical_Variables::penalty * psi(l) * psi(ll) * w * J;
                }
            }
        }
    }

    void getInvCoupledMassPerTime(const DPMVCoupledElement& elem, const Vector<Vector<double>>& couplingMatrix,
                                  const Vector<Vector<double>>& bulkStiffness, Vector<double>& invCMassPerTime)
    {
        // get number of nodes in the bulk element
        const unsigned nnode = bulkStiffness.size();
        // get number of coupled particles in the bulk element
        const unsigned nParticles = invCMassPerTime.size();
        // get coupled particles in the bulk element
        Vector<BaseParticle*> particles;
        if (!useCGMapping())
        { particles = elem.listOfCoupledParticles; }
        else
        { particles = elem.listOfCoupledParticlesExt; }
        // loop over the particles
        for (unsigned m = 0; m < nParticles; m++)
        {
            // get particle mass
            double particleMass = particles[m]->getMass() / particles[m]->getInvCouplingWeight();
            // initialize additional coupling mass from the bulk element
            double cMassPerTimeSquared = 0.0;
            // loop over the nodes
            for (unsigned l = 0; l < nnode; l++)
            {
                // loop over the nodes again
                for (unsigned ll = 0; ll < nnode; ll++)
                {
                    // add contributions to the coupling stiffness or mass over time squared
                    cMassPerTimeSquared += couplingMatrix[l][m] * bulkStiffness[l][ll] * couplingMatrix[ll][m];
                }
            }
            cMassPerTimeSquared *= 0.5;
            // add coupling mass to the particle (in the DPM scale)
            double cMass = cMassPerTimeSquared * Global_Physical_Variables::stiffScale() * pow(getTimeStep(), 2);
            particles[m]->setCouplingMass(cMass);
            // set coupling mass/time only along the diagonal (lumped version)
            double scale = Global_Physical_Variables::massScale() / Global_Physical_Variables::timeScale;
            double oomph_dt = time_pt()->dt(0);
            invCMassPerTime[m] = 1.0 / ( particleMass / getTimeStep() / scale + cMassPerTimeSquared * oomph_dt );
        }
    }

    /*
     * used in computeOneTimeStepForVCoupling to compute weight on particles
     */
    void computeWeightOnParticles()
    {
        // first loop over the coupled bulk elements
        for (DPMVCoupledElement coupled_elem : listOfDPMVCoupledElements_)
        {
            // get pointer to the bulk element
            VELEMENT* elem_pt = coupled_elem.bulk_elem_pt;
            // Find out how many nodes there are
            const unsigned n_node = elem_pt->nnode();
            // Set up memory for the shape/test functions
            Shape psi(n_node);
            // loop over the particles in the bulk element
            for (unsigned i = 0; i < coupled_elem.listOfCoupledParticles.size(); i++)
            {
                auto shapes = coupled_elem.listOfShapesAtParticleCenters[i];
                // Compute weight at particle center
                double weightAtCenter = 0;
                for (unsigned l = 0; l < n_node; l++)
                {
                    CoupledSolidNode* node_pt = dynamic_cast<CoupledSolidNode*>(elem_pt->node_pt(l));
                    weightAtCenter += shapes[l] * node_pt->get_coupling_weight();
                }
                double invWeightAtCenter = 1.0 / weightAtCenter;

                // get point to particle i in the coupled bulk element and set the coupling weight
                auto p = coupled_elem.listOfCoupledParticles[i];
                p->setInvCouplingWeight(invWeightAtCenter);

                // skip particle i if no interactions
                if (p->getInteractions().size() == 0) continue;

                if (explicit_)
                {
                    // loop over interactions with particle i (FIXME accessing coupled interactions is n squared...)
                    for (auto inter : p->getInteractions())
                    {
                        // skip if the coupling weight is already computed
                        if (inter->isCoupled())
                        { continue; }
                        else
                        {
                            // check if a contact point resides within the bulk element
                            Vector<double> s(3, 0.0), x(3, 0.0);
                            GeomObject* geom_obj_pt = 0;
                            // get global coordinates of contact point (DPM)
                            Vec3D xc = inter->getContactPoint();
                            xc /= Global_Physical_Variables::lenScale;
                            x[0] = xc.getX();
                            x[1] = xc.getY();
                            x[2] = xc.getZ();
                            // get local coordinates of contact point (FEM)
                            elem_pt->locate_zeta(x, geom_obj_pt, s);
                            // if true compute weight at contact point
                            if (geom_obj_pt != nullptr)
                            {
                                // Get shape/test functions
                                elem_pt->shape(s, psi);
                                double weightAtContactPoint = 0;
                                for (unsigned l = 0; l < n_node; l++)
                                {
                                    CoupledSolidNode* node_pt = dynamic_cast<CoupledSolidNode*>(elem_pt->node_pt(l));
                                    weightAtContactPoint += node_pt->get_coupling_weight() * psi(l);
                                }
                                inter->setCouplingWeight(weightAtContactPoint);
                                inter->setCoupled();
                            }
                        }
                    }
                }
            }
        }
    }

    void setExplicit(const bool& flag)
    { explicit_ = flag; }

private:
    // vector of bulk elements coupled with discrete particles
    Vector<DPMVCoupledElement> listOfDPMVCoupledElements_;
    // minimal weighting factor
    double weightMin_ = 0.0;
    // flag to use implicit or explicit coupling scheme (VolumeCoupling)
    bool explicit_ = false;
};

#endif //VOLUME_COUPLING_H
