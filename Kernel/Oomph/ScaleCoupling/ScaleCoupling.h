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

#ifndef SCALE_COUPLING_H
#define SCALE_COUPLING_H

#include <utility>

#include "Mercury3D.h"
#include "../Coupling/BaseCoupling.h"
#include "generic.h"
#include "Oomph/RefineableQDPVDElement.h" //
#include "Oomph/SolidProblem.h" //
#include "ScaleCoupledElement.h"
using namespace oomph;

template<class M, class O>
class ScaleCoupling : public BaseCoupling<M,O> {

    /// type of el_pt in O (e.g. RefineableQDPVDElement<3, 2>)
    typedef typename O::ELEMENT ELEMENT;

    /** coupling weight:
     *  0: pure FEM, 1: pure DEM, (0,1): overlap region
     */
    std::function<double(double x,double y,double z)> couplingWeight_;

    /// Stores properties of a coupled element: pointer to the element and list of particles coupled to it.
    struct CoupledElement {
        // pointer to the oomph-lib element
        ELEMENT* el_pt;
        // particles in this element
        std::vector<BaseParticle*> coupledParticles; ///\todo set
    };

    /// Vector of all coupled elements
    std::vector<CoupledElement> coupledElements_;

    /// For all particles, stores coupling properties: coupling force, pointer to coupled element and location in coupled element.
    struct CoupledParticle {
        /// Coupling force
        Vec3D couplingForce;
        /// Element in which this particle resides (null if particle is not coupled)
        CoupledElement* coupledElement_pt = nullptr;
        /// Location of particle in coupled element, returned by locate zeta
        Vector<double> s;

        /// Removes particle from coupledParticle.coupledElement.particles and sets coupledParticle.coupledElement vector to null
        void removeCoupledElement(BaseParticle* particle) {
            if (coupledElement_pt != nullptr) {
                remove(coupledElement_pt->coupledParticles.begin(), coupledElement_pt->coupledParticles.end(), particle);
                coupledElement_pt = nullptr;
                logger(VERBOSE, "Removed particle % from element %", particle->getIndex(), coupledElement_pt);
            }
        }

        /// Remove particle from coupledParticle.coupledElement.particles and set coupledParticle.coupledElement vector to null
        void setCoupledElement(CoupledElement* coupledElementNew, BaseParticle* particle) {
            if (coupledElement_pt != coupledElementNew) {
                if (coupledElement_pt != nullptr) {
                    remove(coupledElement_pt->coupledParticles.begin(), coupledElement_pt->coupledParticles.end(), particle);
                    logger(VERBOSE, "Moved particle % from element % to %", particle->getIndex(), coupledElement_pt, coupledElementNew);
                } else {
                    logger(VERBOSE,"Added particle % to element %", particle->getIndex(), coupledElementNew);
                }
                // add coupledElement to particle
                coupledElement_pt = coupledElementNew;
                // add particle to coupledElement
                coupledElement_pt->coupledParticles.push_back(particle);
            }
        }
    };

    /// The i-th element of this vector describes the coupling properties of the i-th DPM particle
    std::vector<CoupledParticle> coupledParticles_;

    /// Penalty parameter, i.e., proportionality constant between velocity difference and coupling force
    double penalty_ = constants::NaN; ///\todo check if set

public:

    /// get coupled element
    const std::vector<CoupledElement>& getCoupledElements() {
        return coupledElements_;
    }

    /// get coupled particle
    const std::vector<CoupledParticle>& getCoupledParticles() {
        return coupledParticles_;
    }

    /// Sets penalty parameter
    void setPenalty(double penalty) {
        penalty_ = penalty;
    }

    /// Sets weight function (determines which nodes/el_pts are in coupling zone)
    void setCouplingWeight(const std::function<double(double x,double y,double z)>& couplingWeight) {
        couplingWeight_ = couplingWeight;
    }

    /// Computes nStep, the ratio of FEM and DEM time steps, then calls solveSurfaceCoupling(nStep)
    void solveScaleCoupling()
    {
        // compute nStep
        unsigned nStep = O::getOomphTimeStep()/M::getTimeStep();
        if (nStep==0) {
            // if oomph has a smaller time step than Mercury
            logger(INFO, "Set nStep %, change mercuryTimeStep from % to %",
                   nStep, M::getTimeStep(), O::getOomphTimeStep());
            nStep = 1;
            M::setTimeStep(O::getOomphTimeStep());
        } else {
            logger(INFO, "Set nStep %, change oomphTimeStep from % to %",
                   nStep, O::getOomphTimeStep(), nStep * M::getTimeStep());
            O::setOomphTimeStep(nStep * M::getTimeStep());
        }
        // call solve routine
        solveScaleCoupling(nStep);
    }

private:

    /// Solve scale-coupled problem
    void solveScaleCoupling(unsigned nStep)
    {
        // Check main parameters are set
        logger.assert_always(O::getOomphTimeStep()>0, "Oomph time step not initialised");
        logger.assert_always(M::getTimeStep()>0, "Mercury time step not initialised");
        logger.assert_always(std::isfinite(penalty_), "Penalty parameter not set");
        logger.assert_always((bool)couplingWeight_, "Coupling weight is not set");

        // Initialise the coupled element vector
        initialiseCoupledElements();
        
        // This is the main loop advancing over time
        unsigned nDone = 0; //< last written file number
        while (M::getTime() < M::getTimeMax())
        {
            O::actionsBeforeOomphTimeStep();
            // solve the coupled problem for one time step
            computeOneTimeStepScaleCoupling(nStep);
            // write outputs of the oomphProb; this is slaved to the vtk output of Mercury, i.e. an oomph-lib output get written everytime a Mercury vtk file gets written
            if (M::getParticlesWriteVTK() && M::getVtkWriter()->getFileCounter() > nDone)
            {
                O::writeToVTK();
                nDone = M::getVtkWriter()->getFileCounter();
            }
        }

        // Write final output files
        M::finaliseSolve();
    }
    
    /// initialises the coupledElements vector
    void initialiseCoupledElements() {
        // loop through all el_pts, check if any nodes have a fractional coupling weight.
        // If yes, add the el_pt to the coupledElements_
        auto nelement = O::mesh_pt()->nelement();
        for (int e=0; e < nelement; ++e) {
            bool coupledElement = false;
            auto* el_pt = dynamic_cast<ELEMENT*>(O::mesh_pt()->finite_element_pt(e));
            auto nnode = el_pt->nnode();
            Vector<double> nodalCouplingWeights(nnode);
            for (int n=0; n < nnode; ++n) {
                auto* n_pt = dynamic_cast<SolidNode*>(el_pt->node_pt(n));
                nodalCouplingWeights[n] = couplingWeight_(n_pt->xi(0), n_pt->xi(1), n_pt->xi(2));
                if (nodalCouplingWeights[n]>0) { // check for nodalCouplingWeights[n]<1 ?
                    coupledElement = true;
                    logger(VERBOSE,"Coupled element at % % %, weight %",n_pt->xi(0), n_pt->xi(1), n_pt->xi(2), nodalCouplingWeights[n]);
                }
            }
            if (coupledElement) {
                coupledElements_.push_back({el_pt});
                el_pt->set_coupling_weight(nodalCouplingWeights);
            }
        }
        logger(INFO,"% of % el_pts are coupled", coupledElements_.size(), nelement);
    }

    /// What to do in each time step
    void computeOneTimeStepScaleCoupling(unsigned nStep) {
        updateCoupledElements();
        computeCouplingForce();
        this->solveOomph();
        this->solveMercury(nStep);
    };

    /// Updates the coupledElements and coupledParticles vectors: finds which particles interact with which element
    void updateCoupledElements()
    {
        /// 1) Find min and max of coupled region
        Vector<double> min {constants::inf, constants::inf, constants::inf};
        Vector<double> max {-constants::inf, -constants::inf, -constants::inf};
        // loop over all coupled elements
        for (auto& coupledElement : coupledElements_)
        {
            ELEMENT* el_pt = coupledElement.el_pt;
            const unsigned nnode = el_pt->nnode();
            const unsigned dim = el_pt->dim(); // /todo I need to use el_pt->ndim=3, not el_pt->ndim=0, why?
            // Find min and max
            for (int n=0; n<nnode; ++n) {
                auto* n_pt = dynamic_cast<SolidNode *>(el_pt->node_pt(n));
                for (int d = 0; d < dim; ++d) {
                    min[d] = std::min(min[d], n_pt->x(d));
                    max[d] = std::max(max[d], n_pt->x(d));
                }
            }
        }
        logger(VERBOSE,"Coupling zone: % < x < %, % < y < %, % < z < % ",min[0],max[0],min[1],max[1],min[2],max[2]);

        /// 2) update the coupledParticles vector, and the particles vector of each coupledElement
        // resize the coupledParticles vector
        coupledParticles_.resize(M::particleHandler.getSize());
        // now set the coupledElement pointer for each coupledParticles, and the particles vector for each coupledElement
        for (int p=0; p<coupledParticles_.size(); ++p) {
            CoupledParticle& coupledParticle = coupledParticles_[p];
            Vec3D pos = M::particleHandler.getObject(p)->getPosition();
            // this is the pointer we need to set
            if (pos.X>=min[0] && pos.X<=max[0] && pos.Y>=min[1] && pos.Y<=max[1] && pos.Z>=min[2] && pos.Z<=max[2]) {
                logger(VERBOSE,"Particle % of % is in coupling zone",p, M::particleHandler.getSize());
                // if near the coupling zone
                Vector<double> x {pos.X, pos.Y, pos.Z};
                Vector<double>& s = coupledParticle.s;
                GeomObject* geom_obj_pt;
                // if el_pt is already set and still valid, keep
                if (coupledParticle.coupledElement_pt != nullptr) {
                    coupledParticle.coupledElement_pt->el_pt->locate_zeta(x, geom_obj_pt, s);
                    if (geom_obj_pt!=nullptr) {
                        logger(VERBOSE,"Particle % remains coupled with element %",p, coupledParticle.coupledElement_pt->el_pt);
                        continue;
                    }
                }
                // otherwise, search for the right element 
                for (auto& coupledElement : coupledElements_) {
                    coupledElement.el_pt->locate_zeta(x,geom_obj_pt,s);
                    if (geom_obj_pt!=nullptr) {
                        logger(VERBOSE,"Coupled particle % to element % at % % %",p, &coupledElement, s[0], s[1], s[2]);
                        //\todo If a containing element is found, we stop; note though, particles on the boundary could belong to multiple elements; also not clear to me whether force should be applied to global basis function or local
                        coupledParticle.setCoupledElement(&coupledElement, M::particleHandler.getObject(p));
                        break;
                    }
                }
            } else {
                // if away from the coupling, set coupledElement to null
                if (coupledParticle.coupledElement_pt) logger(VERBOSE, "Removed particle % from element %", p, coupledParticle.coupledElement_pt);
                coupledParticle.removeCoupledElement(M::particleHandler.getObject(p));
            }
        }
    }

    /// Computes projectionMatrix
    Vector<Vector<double>> computeProjectionMatrix(const CoupledElement& coupledElement)
    {
        Vector<Vector<double>> projectionMatrix;
        // get constants
        const unsigned nParticles = coupledElement.coupledParticles.size();
        const unsigned nnode = coupledElement.el_pt->nnode();
        const unsigned dim = coupledElement.el_pt->dim();

        // prepare projection matrix
        projectionMatrix.resize(nnode,Vector<double>(nParticles,0.0));

        /// Construct mapping with FEM basis functions
        // storage for position s and shape psi
        Vector<double> pos(dim,0.0);
        Shape psi(nnode);
        // loop over shape functions at the nodes
        // loop over shape functions at the particles
        for (unsigned p = 0; p < nParticles; p++) {
            BaseParticle* particle = coupledElement.coupledParticles[p];
            // get the shape function evaluated at the center of particle p
            coupledElement.el_pt->shape(coupledParticles_[particle->getIndex()].s, psi);
            for (unsigned n=0; n < nnode; n++) {
                projectionMatrix[n][p] = psi(n) * particle->getVolume();
            }
        }
        // normalize each row of the projection matrix to sum to one
        for (unsigned n=0; n < nnode; n++) {
            double sum = 0;
            for (unsigned p = 0; p<nParticles; p++) {
                sum += projectionMatrix[n][p];
            }
            if (sum != 0) {
                for (unsigned p = 0; p < nParticles; p++) {
                    projectionMatrix[n][p] /= sum;
                }
            }
        }
        return projectionMatrix;
    }

    /// Computes coupling force for each element and particle
    void computeCouplingForce()
    {
        // first loop over the coupled bulk el_pts
        for (const auto& coupledElement : coupledElements_)
        {
            // construct projection matrix, defined as projectionMatrix_{n,p} = N_{n,p}*V_p / sum_p(N_{n,p}*V_p) 
            // (shape function N, volume V, node n, particle p)
            Vector<Vector<double>> projectionMatrix = computeProjectionMatrix(coupledElement);
            // get pointer to the bulk el_pt
            ELEMENT* el_pt = coupledElement.el_pt;
            const unsigned nParticles = coupledElement.coupledParticles.size();
            // get number of nodes in the bulk el_pt
            const unsigned nnode = coupledElement.el_pt->nnode();
            // prepare vector of nodal velocities projected from particle centers
            const unsigned dim = el_pt->dim();
            // get coupling force at nodal positions, based on penalizing the difference in velocity
            computeNodalCouplingForces(coupledElement, projectionMatrix);
            // project coupling forces from nodes to particles (note projectionMatrix^{-1} = projectionMatrix^T)
            for (unsigned p=0; p < nParticles; p++) {
                BaseParticle* particle = coupledElement.coupledParticles[p];
                CoupledParticle& coupledParticle = coupledParticles_[particle->getIndex()];
                for (unsigned d=0; d < dim; d++) {
                    double forceComponent = 0;
                    for (unsigned n=0; n < nnode; n++) {
                        // note projectionMatrix^{-1} = projectionMatrix^T
                        forceComponent += projectionMatrix[n][p] * el_pt->get_coupling_residual(n, d);
                    }
                    coupledParticle.couplingForce.setComponent(d,forceComponent);
                }
                logger(VERBOSE,"Coupling force % on particle %", coupledParticle.couplingForce,p);
            }
        }
    }

    /// Computes coupling force for each node
    void computeNodalCouplingForces(const CoupledElement& coupledElement, const Vector<Vector<double>>& projectionMatrix)
    {
        // return if no particle in cell
        if (coupledElement.coupledParticles.empty()) return;

        // a few shortcut variables
        const ELEMENT* el_pt = coupledElement.el_pt;
        const unsigned nnode = el_pt->nnode();
        const unsigned dim = el_pt->dim();
        const unsigned nParticles = coupledElement.coupledParticles.size();
        
        /// 1) compute velocity difference at nodal positions
        Vector<Vector<double>> nodalVelocityDifference(nnode,Vector<double>(dim,0.0));
        // get projected nodal velocity field from particle velocities
        for (unsigned n=0; n < nnode; n++) {
            for (unsigned d=0; d < dim; d++) {
                double nodalVelocityDPM = 0;
                for (unsigned p=0; p < nParticles; p++) {
                    nodalVelocityDPM += projectionMatrix[n][p] * coupledElement.coupledParticles[p]->getVelocity().getComponent(d);
                }
                double nodalVelocityFEM = el_pt->dnodal_position_gen_dt(n,0,d);
                nodalVelocityDifference[n][d] = nodalVelocityDPM-nodalVelocityFEM;
            }
        }

        /// 2) compute nodal coupling force, penalizing the difference in velocity
        Vector<Vector<double> > nodalCouplingForces(nnode, Vector<double>(dim, 0.0));
        
        // Set up memory for the shape functions
        Shape psi(nnode);
        DShape dpsidxi(nnode,dim);
        const unsigned n_in = el_pt->integral_pt()->nweight();
        // Loop over the integration points
        for (unsigned in=0; in<n_in; in++)
        {
            // Get the integral weight
            double w = el_pt->integral_pt()->weight(in);
            // Call the derivatives of the shape functions (and get Jacobian)
            double J = el_pt->dshape_lagrangian_at_knot(in,psi,dpsidxi);
            // Loop over the test functions, nodes of the element
            for (unsigned n=0; n < nnode; n++)
            {
                // Loop over the nodes of the element again
                for (unsigned nn=0; nn < nnode; nn++)
                {
                    double volume = penalty_ * psi(nn) * psi(n) * w * J;

                    for (unsigned d=0; d < dim; d++)
                    {
                        //nodalVelocityDifference[n][d] = 1/penalty_/time_pt()->dt(0);
                        double displacement = nodalVelocityDifference[n][d] * O::time_pt()->dt(0);
                        nodalCouplingForces[n][d] -= volume * displacement;
                        //logger(INFO, "w % J % psi % % displacement % nodalCouplingForces % dt %", w, J, psi(nn), psi(n), displacement, nodalCouplingForces[n][d], time_pt()->dt(0));
                    }
                }
            }
        }

        /// 3) assign it to the coupled bulk element to be used by ELEMENT::fill_in_contribution_to_residuals(...)
        coupledElement.el_pt->set_coupling_residual(nodalCouplingForces);

        logger(VERBOSE,"VelocityDifference % % %, coupling force % % % on element %, node 0, particles %",
               nodalVelocityDifference[0][0], nodalVelocityDifference[0][1], nodalVelocityDifference[0][2],
               nodalCouplingForces[0][0], nodalCouplingForces[0][1], nodalCouplingForces[0][2],
               el_pt, coupledElement.coupledParticles.size());
    }

    /// Applies coupling force to MercuryDPM in each time step
    void computeExternalForces(BaseParticle* p) override
    {
        if (!p->isFixed())
        {
            // Applying the force due to gravity (F = m.g)
            p->addForce(M::getGravity() * p->getMass());
            // add particle coupling force 1/w_i*f_i^c to the total force
            const CoupledParticle& cp = coupledParticles_[p->getIndex()];
            double couplingWeight = couplingWeight_(p->getPosition().X,p->getPosition().Y,p->getPosition().Z);
            p->addForce(cp.couplingForce / couplingWeight);
        }
    }
};

#endif //SCALE_COUPLING_H
