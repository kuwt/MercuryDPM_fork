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

#ifndef SCALECOUPLEDELEMENT_H
#define SCALECOUPLEDELEMENT_H
#include "Logger.h"
#include "generic.h"
//#include "Oomph/Coupling/CoupledSolidNodes.h"

namespace oomph
{

/// Wrapper class for solid elements to be coupled with discrete solid particles on the surfaces
template<class ELEMENT>
class ScaleCoupledElement : public ELEMENT {
public:

    /// Empties the coupling residual
    void clear_coupling_residual() {
        coupling_residual.clear();
    }

    /// Sets the coupling residual
    void set_coupling_residual(Vector<Vector<double>> &residual) {
        coupling_residual = residual;
    }

    /// Returns the coupling residual
    double get_coupling_residual(const unsigned& n, const unsigned& d)
    {
        return coupling_residual.empty()?0.0:coupling_residual[n][d];
    }

    /// Empties the coupling Jacobian
    void clear_coupling_jacobian() {
        coupling_jacobian.clear();
    }

    /// Sets the coupling residual
    void set_coupling_jacobian(Vector<Vector<double>> &jacobian) {
        coupling_jacobian = jacobian;
    }

    /// Sets the coupling weight
    void set_coupling_weight(Vector<double> &weight) {
        coupling_weight = weight;
    }

    /// Sets the coupling weight
    double get_coupling_weight(const unsigned& n) {
        return coupling_weight.empty()?0.0:coupling_weight[n];
    }

    /// Adds the elements contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double> &residuals) override
    {
        // Call the generic residuals function with flag set to 0 using a dummy matrix argument
        fill_in_generic_contribution_to_residuals_pvd(residuals,GeneralisedElement::Dummy_matrix,0);

        // Add nodal coupling force to the residuals
        add_internal_coupling_forces_to_residuals(residuals,GeneralisedElement::Dummy_matrix,0);

        //if (!coupling_residual.empty()) logger(INFO,"fill_in_contribution_to_residuals, element %", this);
    }

    /// Add the elements contribution to its residual vector and Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double> &residuals, DenseMatrix<double> &jacobian) override
    {
        //Call the generic routine with the flag set to 1
        fill_in_generic_contribution_to_residuals_pvd(residuals,jacobian,1);

        // Add nodal coupling force to the residuals
        add_internal_coupling_forces_to_residuals(residuals,jacobian,1);

        //if (!coupling_residual.empty()) logger(INFO,"fill_in_contribution_to_jacobian, element %", this);
    }

    /// \todo TW I think this just addes the coupling weight
    void fill_in_generic_contribution_to_residuals_pvd(Vector<double> &residuals,
                                                       DenseMatrix<double> &jacobian,
                                                       const unsigned& flag) override
    {

        // Get problem dimension
        const unsigned DIM = this->dim();

        // Simply set up initial condition?
        if (this->Solid_ic_pt!=0)
        {
            this->fill_in_residuals_for_solid_ic(residuals);
            return;
        }

        //Find out how many nodes there are
        const unsigned n_node = this->nnode();

        //Find out how many positional dofs there are
        const unsigned n_position_type = this->nnodal_position_type();

        //Set up memory for the shape functions
        Shape psi(n_node,n_position_type);
        DShape dpsidxi(n_node,n_position_type,DIM);

        //Set the value of Nintpt -- the number of integration points
        const unsigned n_intpt = this->integral_pt()->nweight();

        //Set the vector to hold the local coordinates in the element
        Vector<double> s(DIM);

        // Timescale ratio (non-dim density)
        double lambda_sq = this->lambda_sq();

        // Time factor
        double time_factor=0.0;
        if (lambda_sq>0)
        {
            time_factor=this->node_pt(0)->position_time_stepper_pt()->weight(2,0);
        }

        //Integer to store the local equation number
        int local_eqn=0;

        //Loop over the integration points
        for(unsigned ipt=0;ipt<n_intpt;ipt++)
        {
            //Assign the values of s
            for(unsigned i=0;i<DIM;++i) {s[i] = this->integral_pt()->knot(ipt,i);}

            //Get the integral weight
            double w = this->integral_pt()->weight(ipt);

            //Call the derivatives of the shape functions (and get Jacobian)
            double J = this->dshape_lagrangian_at_knot(ipt,psi,dpsidxi);

            //Get the coupling weight at the integration point
            double coupling_w = 0;
            if (!coupling_weight.empty()) {
                for (unsigned l = 0; l < n_node; l++) {
                    double psi_ = psi(l);
                    coupling_w += psi_ * coupling_weight[l];
                }
            }

            //Calculate interpolated values of the derivative of global position
            //wrt lagrangian coordinates
            DenseMatrix<double> interpolated_G(DIM);

            // Setup memory for accelerations
            Vector<double> accel(DIM);

            //Initialise to zero
            for(unsigned i=0;i<DIM;i++)
            {
                // Initialise acclerations
                accel[i]=0.0;
                for(unsigned j=0;j<DIM;j++)
                {
                    interpolated_G(i,j) = 0.0;
                }
            }

            //Storage for Lagrangian coordinates (initialised to zero)
            Vector<double> interpolated_xi(DIM,0.0);

            //Calculate displacements and derivatives and lagrangian coordinates
            for(unsigned l=0;l<n_node;l++)
            {
                //Loop over positional dofs
                for(unsigned k=0;k<n_position_type;k++)
                {
                    double psi_ = psi(l,k);
                    //Loop over displacement components (deformed position)
                    for(unsigned i=0;i<DIM;i++)
                    {
                        //Calculate the Lagrangian coordinates and the accelerations
                        interpolated_xi[i] += this->lagrangian_position_gen(l,k,i)*psi_;

                        // Only compute accelerations if inertia is switched on
                        if ((lambda_sq>0.0)&&(this->Unsteady))
                        {
                            accel[i] += this->dnodal_position_gen_dt(2,l,k,i)*psi_;
                        }

                        //Loop over derivative directions
                        for(unsigned j=0;j<DIM;j++)
                        {
                            interpolated_G(j,i) +=
                                this->nodal_position_gen(l,k,i)*dpsidxi(l,k,j);
                        }
                    }
                }
            }

            //Get isotropic growth factor
            double gamma=1.0;
            this->get_isotropic_growth(ipt,s,interpolated_xi,gamma);


            //Get body force at current time
            Vector<double> b(DIM);
            this->body_force(interpolated_xi,b);

            // We use Cartesian coordinates as the reference coordinate
            // system. In this case the undeformed metric tensor is always
            // the identity matrix -- stretched by the isotropic growth
            double diag_entry=pow(gamma,2.0/double(DIM));
            DenseMatrix<double> g(DIM);
            for(unsigned i=0;i<DIM;i++)
            {
                for(unsigned j=0;j<DIM;j++)
                {
                    if(i==j) {g(i,j) = diag_entry;}
                    else {g(i,j) = 0.0;}
                }
            }

            //Premultiply the undeformed volume ratio (from the isotropic
            // growth), the integral weights, the coupling weights, and the Jacobian
            double W = gamma*w*(1.0-coupling_w)*J;

            //Declare and calculate the deformed metric tensor
            DenseMatrix<double> G(DIM);

            //Assign values of G
            for(unsigned i=0;i<DIM;i++)
            {
                //Do upper half of matrix
                for(unsigned j=i;j<DIM;j++)
                {
                    //Initialise G(i,j) to zero
                    G(i,j) = 0.0;
                    //Now calculate the dot product
                    for(unsigned k=0;k<DIM;k++)
                    {
                        G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
                    }
                }
                //Matrix is symmetric so just copy lower half
                for(unsigned j=0;j<i;j++)
                {
                    G(i,j) = G(j,i);
                }
            }

            //Now calculate the stress tensor from the constitutive law
            DenseMatrix<double> sigma(DIM);
            ELEMENT::get_stress(g,G,sigma);

            // Add pre-stress
            for(unsigned i=0;i<DIM;i++)
            {
                for(unsigned j=0;j<DIM;j++)
                {
                    sigma(i,j) += this->prestress(i,j,interpolated_xi);
                }
            }

            // Get stress derivative by FD only needed for Jacobian
            //-----------------------------------------------------

            // Stress derivative
            RankFourTensor<double> d_stress_dG(DIM,DIM,DIM,DIM,0.0);
            // Derivative of metric tensor w.r.t. to nodal coords
            RankFiveTensor<double> d_G_dX(n_node,n_position_type,DIM,DIM,DIM,0.0);

            // Get Jacobian too?
            if (flag==1)
            {
                // Derivative of metric tensor w.r.t. to discrete positional dofs
                // NOTE: Since G is symmetric we only compute the upper triangle
                //          and DO NOT copy the entries across. Subsequent computations
                //          must (and, in fact, do) therefore only operate with upper
                //          triangular entries
                for (unsigned ll=0;ll<n_node;ll++)
                {
                    for (unsigned kk=0;kk<n_position_type;kk++)
                    {
                        for (unsigned ii=0;ii<DIM;ii++)
                        {
                            for (unsigned aa=0;aa<DIM;aa++)
                            {
                                for (unsigned bb=aa;bb<DIM;bb++)
                                {
                                    d_G_dX(ll,kk,ii,aa,bb)=
                                        interpolated_G(aa,ii)*dpsidxi(ll,kk,bb)+
                                        interpolated_G(bb,ii)*dpsidxi(ll,kk,aa);
                                }
                            }
                        }
                    }
                }

                //Get the "upper triangular" entries of the derivatives of the stress
                //tensor with respect to G
                this->get_d_stress_dG_upper(g,G,sigma,d_stress_dG);
            }

//=====EQUATIONS OF ELASTICITY FROM PRINCIPLE OF VIRTUAL DISPLACEMENTS========

            //Loop over the test functions, nodes of the element
            for(unsigned l=0;l<n_node;l++)
            {
                //Loop of types of dofs
                for(unsigned k=0;k<n_position_type;k++)
                {
                    // Offset for faster access
                    const unsigned offset5=dpsidxi.offset(l ,k);

                    //Loop over the displacement components
                    for(unsigned i=0;i<DIM;i++)
                    {
                        //Get the equation number
                        local_eqn = this->position_local_eqn(l,k,i);

                        /*IF it's not a boundary condition*/
                        if(local_eqn >= 0)
                        {
                            //Initialise contribution to sum
                            double sum=0.0;

                            // Acceleration and body force
                            sum+=(lambda_sq*accel[i]-b[i])*psi(l,k);

                            // Stress term
                            for(unsigned a=0;a<DIM;a++)
                            {
                                unsigned count=offset5;
                                for(unsigned b=0;b<DIM;b++)
                                {
                                    //Add the stress terms to the residuals
                                    sum+=sigma(a,b)*interpolated_G(a,i)*
                                         dpsidxi.raw_direct_access(count);
                                    ++count;
                                }
                            }
                            residuals[local_eqn] += W*sum;

                            // Get Jacobian too?
                            if (flag==1)
                            {
                                // Offset for faster access in general stress loop
                                const unsigned offset1=d_G_dX.offset( l, k, i);

                                //Loop over the nodes of the element again
                                for(unsigned ll=0;ll<n_node;ll++)
                                {
                                    //Loop of types of dofs again
                                    for(unsigned kk=0;kk<n_position_type;kk++)
                                    {
                                        //Loop over the displacement components again
                                        for(unsigned ii=0;ii<DIM;ii++)
                                        {
                                            //Get the number of the unknown
                                            int local_unknown = this->position_local_eqn(ll,kk,ii);

                                            /*IF it's not a boundary condition*/
                                            if(local_unknown >= 0)
                                            {
                                                // Offset for faster access in general stress loop
                                                const unsigned offset2=d_G_dX.offset( ll, kk, ii);
                                                const unsigned offset4=dpsidxi.offset(ll, kk);

                                                // General stress term
                                                //--------------------
                                                double sum=0.0;
                                                unsigned count1=offset1;
                                                for(unsigned a=0;a<DIM;a++)
                                                {
                                                    // Bump up direct access because we're only
                                                    // accessing upper triangle
                                                    count1+=a;
                                                    for(unsigned b=a;b<DIM;b++)
                                                    {
                                                        double factor=d_G_dX.raw_direct_access(count1);
                                                        if (a==b) factor*=0.5;

                                                        // Offset for faster access
                                                        unsigned offset3=d_stress_dG.offset(a,b);
                                                        unsigned count2=offset2;
                                                        unsigned count3=offset3;

                                                        for(unsigned aa=0;aa<DIM;aa++)
                                                        {
                                                            // Bump up direct access because we're only
                                                            // accessing upper triangle
                                                            count2+=aa;
                                                            count3+=aa;

                                                            // Only upper half of derivatives w.r.t. symm tensor
                                                            for(unsigned bb=aa;bb<DIM;bb++)
                                                            {
                                                                sum+=factor*
                                                                     d_stress_dG.raw_direct_access(count3)*
                                                                     d_G_dX.raw_direct_access(count2);
                                                                ++count2;
                                                                ++count3;
                                                            }
                                                        }
                                                        ++count1;
                                                    }

                                                }

                                                // Multiply by weight and add contribution
                                                // (Add directly because this bit is nonsymmetric)
                                                jacobian(local_eqn,local_unknown)+=sum*W;

                                                // Only upper triangle (no separate test for bc as
                                                // local_eqn is already nonnegative)
                                                if((i==ii) &&   (local_unknown >= local_eqn))
                                                {
                                                    //Initialise contribution
                                                    double sum=0.0;

                                                    // Inertia term
                                                    sum+=lambda_sq*time_factor*psi(ll,kk)*psi(l,k);

                                                    // Stress term
                                                    unsigned count4=offset4;
                                                    for(unsigned a=0;a<DIM;a++)
                                                    {
                                                        //Cache term
                                                        const double factor=
                                                            dpsidxi.raw_direct_access(count4);// ll ,kk
                                                        ++count4;

                                                        unsigned count5=offset5;
                                                        for(unsigned b=0;b<DIM;b++)
                                                        {
                                                            sum+=sigma(a,b)*factor*
                                                                 dpsidxi.raw_direct_access(count5); // l   ,k
                                                            ++count5;
                                                        }
                                                    }
                                                    //Add contribution to jacobian
                                                    jacobian(local_eqn,local_unknown) += sum*W;
                                                    //Add to lower triangular section
                                                    if(local_eqn != local_unknown)
                                                    {
                                                        jacobian(local_unknown,local_eqn) += sum*W;
                                                    }
                                                }

                                            } //End of if not boundary condition
                                        }
                                    }
                                }
                            }

                        } //End of if not boundary condition

                    } //End of loop over coordinate directions
                } //End of loop over type of dof
            } //End of loop over shape functions
        } //End of loop over integration points
    }

    /// Add the point source contribution to the residual vector
    void add_internal_coupling_forces_to_residuals(Vector<double> &residuals,
                                                   DenseMatrix<double> &jacobian,
                                                   const unsigned& flag)
    {
        // No further action if no coupling forces
        if (coupling_residual.empty()) return;

        // Find out how many nodes there are
        const unsigned n_node = this->nnode();
        const unsigned DIM = this->dim();

        // Find out how many positional dofs there are
        unsigned n_position_type = this->nnodal_position_type();

        // Set up memory for the shape/test functions
        Shape psi(n_node);

        // Loop over all nodes belonging to the element
        for(unsigned l=0;l<n_node;l++)
        {
            // Loop of types of dofs
            for(unsigned k=0; k<n_position_type; k++)
            {
                // Loop over the force components
                for(unsigned i=0; i<DIM; i++)
                {
                    // Get the local equation
                    int local_eqn = this->position_local_eqn(l,k,i);
                    // IF it's not a boundary condition and the interpolated force is nonzero
                    if (local_eqn >= 0)
                    {
                        // add the nodal coupling residual to the global residual vector
                        residuals[local_eqn] += coupling_residual[l][i];
                    }
                    // Get Jacobian too?
                    if (flag==1)
                    {
                        // No further action if no coupling forces
                        if (coupling_jacobian.empty()) continue;

                        // Loop over the nodes of the element again
                        for(unsigned ll=0;ll<n_node;ll++)
                        {
                            // Loop of types of dofs again
                            for(unsigned kk=0;kk<n_position_type;kk++)
                            {
                                // Loop over the force components again
                                for(unsigned ii=0;ii<DIM;ii++)
                                {
                                    // Get the number of the unknown
                                    int local_unknown = this->position_local_eqn(ll,kk,ii);
                                    // IF it's not a boundary condition
                                    if ((local_unknown >= 0) && (i==ii) && (local_unknown >= local_eqn))
                                    {
                                        // add the nodal coupling jacobian to the global jacobian matrix
                                        jacobian(local_eqn,local_unknown) += coupling_jacobian[l][ll];
                                        if(local_eqn != local_unknown)
                                        {
                                            jacobian(local_unknown,local_eqn) += coupling_jacobian[l][ll];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

private:
    /// Contribution of coupling to residual (residuals[local_eqn] += coupling_residual[node][dim];)
    // For comparison residuals[local_eqn] +=gamma*w*(1.0-coupling_w)*J*b[i]*psi(nn,d);
    // penalty_ * psi(nn) * psi(n) * w * J * nodalVelocityDifference[n][d] * time_pt()->dt(0)
    Vector<Vector<double>> coupling_residual;
    // Force on particle:
    // force[d] += projectionMatrix[n][p] * el_pt->get_coupling_residual(n, d);

    /// Contribution of coupling to jacobian (jacobian(local_eqn,local_unknown) += coupling_jacobian[node0][node1];)
    Vector<Vector<double>> coupling_jacobian;

    /// Nodal coupling weight (interpolated to compute coupling weight at integration points)
    Vector<double> coupling_weight;
};


/// FaceGeometry of wrapped element is the same as the underlying element
template<class ELEMENT>
class FaceGeometry<ScaleCoupledElement<ELEMENT>> :
    public virtual FaceGeometry<ELEMENT> {
public:
    //FaceGeometry() : FaceGeometry<ELEMENT>() {}
};

}

#endif //SCALECOUPLEDELEMENT_H
