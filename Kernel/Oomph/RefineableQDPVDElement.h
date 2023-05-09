//
// Created by Thomas Weinhart on 29/04/2022.
//

#ifndef MERCURYDPM_REFINEABLEQDPVDELEMENT_H
#define MERCURYDPM_REFINEABLEQDPVDELEMENT_H
#include "refineable_solid_elements.h"

namespace oomph
{

/**
 * Adds a dissipative force to the solid constitutive equation
 * @tparam DIM spatial dimension of the problem
 * @tparam NNODE_1D order of fem space +1
 */
template<unsigned DIM, unsigned NNODE_1D>
class RefineableQDPVDElement : public RefineableQPVDElement<DIM,NNODE_1D>
{
public:

    RefineableQDPVDElement() : RefineableQPVDElement<DIM,NNODE_1D> () {}

    void fill_in_generic_contribution_to_residuals_pvd(
        Vector<double>& residuals,
        DenseMatrix<double>& jacobian,
        const unsigned& flag )

    {
#ifdef PARANOID
        // Check if the constitutive equation requires the explicit imposition of an
        // incompressibility constraint
        if (this->Constitutive_law_pt->requires_incompressibility_constraint())
        {
            throw OomphLibError("RefineablePVDEquations cannot be used with "
                                 "incompressible constitutive laws.",
                                 OOMPH_CURRENT_FUNCTION,
                                 OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Simply set up initial condition?
        if (this->Solid_ic_pt != 0)
        {
            this->get_residuals_for_solid_ic(residuals);
            return;
        }

        // Find out how many nodes there are
        const unsigned n_node = this->nnode();

        // Find out how many positional dofs there are
        const unsigned n_position_type = this->nnodal_position_type();

        // Integers to store local equation numbers
        int local_eqn = 0;

        // Timescale ratio (non-dim density)
        double lambda_sq = this->lambda_sq();

        // Time factor
        double time_factor = 0.0;
        double time_factor1 = 0.0;
        if(lambda_sq>0)
        {
            time_factor = this->node_pt(0)->position_time_stepper_pt()->weight(2, 0);
            time_factor1 = this->node_pt(0)->position_time_stepper_pt()->weight(1, 0);
        }

        // Set up memory for the shape functions
        Shape psi(n_node, n_position_type);
        DShape dpsidxi(n_node, n_position_type, DIM);

        // Set the value of n_intpt -- the number of integration points
        const unsigned n_intpt = this->integral_pt()->nweight();

        // Set the vector to hold the local coordinates in the element
        Vector<double> s(DIM);

        // Loop over the integration points
        for (unsigned ipt = 0; ipt < n_intpt; ipt++)
        {
            // Assign the values of s
            for (unsigned i = 0; i < DIM; ++i)
            {
                s[i] = this->integral_pt()->knot(ipt, i);
            }

            // Get the integral weight
            double w = this->integral_pt()->weight(ipt);

            // Call the derivatives of the shape functions (and get Jacobian)
            double J = this->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);

            // Calculate interpolated values of the derivative of global position
            // wrt lagrangian coordinates
            DenseMatrix<double> interpolated_G(DIM, DIM, 0.0);

            // Setup memory for accelerations
            Vector<double> accel(DIM, 0.0);

            // Storage for Lagrangian coordinates (initialised to zero)
            Vector<double> interpolated_xi(DIM, 0.0);

            // Calculate displacements and derivatives and lagrangian coordinates
            for (unsigned l = 0; l < n_node; l++)
            {
                // Loop over positional dofs
                for (unsigned k = 0; k < n_position_type; k++)
                {
                    double psi_ = psi(l, k);
                    // Loop over displacement components (deformed position)
                    for (unsigned i = 0; i < DIM; i++)
                    {
                        // Calculate the Lagrangian coordinates and the accelerations
                        interpolated_xi[i] += this->lagrangian_position_gen(l, k, i) * psi_;

                        // Only compute accelerations if inertia is switched on
                        // otherwise the timestepper might not be able to
                        // work out dx_gen_dt(2,...)
                        if ((lambda_sq > 0.0) && (this->Unsteady))
                        {
                            accel[i] += this->dnodal_position_gen_dt(2, l, k, i) * psi_;
                            accel[i] += this->dissipation_ *
                                this->dnodal_position_gen_dt(1, l, k, i) * psi_; //TW
                        }

                        // Loop over derivative directions
                        for (unsigned j = 0; j < DIM; j++)
                        {
                            interpolated_G(j, i) +=
                                this->nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
                        }
                    }
                }
            }

            // Get isotropic growth factor
            double gamma = 1.0;
            this->get_isotropic_growth(ipt, s, interpolated_xi, gamma);

            // Get body force at current time
            Vector<double> b(DIM);
            this->body_force(interpolated_xi, b);

            // We use Cartesian coordinates as the reference coordinate
            // system. In this case the undeformed metric tensor is always
            // the identity matrix -- stretched by the isotropic growth
            double diag_entry = pow(gamma, 2.0 / double(DIM));
            DenseMatrix<double> g(DIM);
            for (unsigned i = 0; i < DIM; i++)
            {
                for (unsigned j = 0; j < DIM; j++)
                {
                    if (i == j)
                    {
                        g(i, j) = diag_entry;
                    }
                    else
                    {
                        g(i, j) = 0.0;
                    }
                }
            }

            // Premultiply the undeformed volume ratio (from the isotropic
            // growth), the weights and the Jacobian
            double W = gamma * w * J;

            // Declare and calculate the deformed metric tensor
            DenseMatrix<double> G(DIM);

            // Assign values of G
            for (unsigned i = 0; i < DIM; i++)
            {
                // Do upper half of matrix
                for (unsigned j = i; j < DIM; j++)
                {
                    // Initialise G(i,j) to zero
                    G(i, j) = 0.0;
                    // Now calculate the dot product
                    for (unsigned k = 0; k < DIM; k++)
                    {
                        G(i, j) += interpolated_G(i, k) * interpolated_G(j, k);
                    }
                }
                // Matrix is symmetric so just copy lower half
                for (unsigned j = 0; j < i; j++)
                {
                    G(i, j) = G(j, i);
                }
            }

            // Now calculate the stress tensor from the constitutive law
            DenseMatrix<double> sigma(DIM);
            this->get_stress(g, G, sigma);

            // Get stress derivative by FD only needed for Jacobian
            //-----------------------------------------------------

            // Stress derivative
            RankFourTensor<double> d_stress_dG(DIM, DIM, DIM, DIM, 0.0);
            // Derivative of metric tensor w.r.t. to nodal coords
            RankFiveTensor<double> d_G_dX(
                n_node, n_position_type, DIM, DIM, DIM, 0.0);

            // Get Jacobian too?
            if (flag == 1)
            {
                // Derivative of metric tensor w.r.t. to discrete positional dofs
                // NOTE: Since G is symmetric we only compute the upper triangle
                //       and DO NOT copy the entries across. Subsequent computations
                //       must (and, in fact, do) therefore only operate with upper
                //       triangular entries
                for (unsigned ll = 0; ll < n_node; ll++)
                {
                    for (unsigned kk = 0; kk < n_position_type; kk++)
                    {
                        for (unsigned ii = 0; ii < DIM; ii++)
                        {
                            for (unsigned aa = 0; aa < DIM; aa++)
                            {
                                for (unsigned bb = aa; bb < DIM; bb++)
                                {
                                    d_G_dX(ll, kk, ii, aa, bb) =
                                        interpolated_G(aa, ii) * dpsidxi(ll, kk, bb) +
                                        interpolated_G(bb, ii) * dpsidxi(ll, kk, aa);
                                }
                            }
                        }
                    }
                }

                // Get the "upper triangular"
                // entries of the derivatives of the stress tensor with
                // respect to G
                this->get_d_stress_dG_upper(g, G, sigma, d_stress_dG);
            }


            // Add pre-stress
            for (unsigned i = 0; i < DIM; i++)
            {
                for (unsigned j = 0; j < DIM; j++)
                {
                    sigma(i, j) += this->prestress(i, j, interpolated_xi);
                }
            }

            //=====EQUATIONS OF ELASTICITY FROM PRINCIPLE OF VIRTUAL
            // DISPLACEMENTS========


            // Default setting for non-hanging node
            unsigned n_master = 1;
            double hang_weight = 1.0;

            // Loop over the test functions, nodes of the element
            for (unsigned l = 0; l < n_node; l++)
            {
                // Get pointer to local node l
                Node* local_node_pt = this->node_pt(l);

                // Cache hang status
                bool is_hanging = local_node_pt->is_hanging();

                // If the node is a hanging node
                if (is_hanging)
                {
                    n_master = local_node_pt->hanging_pt()->nmaster();
                }
                // Otherwise the node is its own master
                else
                {
                    n_master = 1;
                }


                // Storage for local equation numbers at node indexed by
                // type and direction
                DenseMatrix<int> position_local_eqn_at_node(n_position_type, DIM);

                // Loop over the master nodes
                for (unsigned m = 0; m < n_master; m++)
                {
                    if (is_hanging)
                    {
                        // Find the equation numbers
                        position_local_eqn_at_node = this->local_position_hang_eqn(
                            local_node_pt->hanging_pt()->master_node_pt(m));

                        // Find the hanging node weight
                        hang_weight = local_node_pt->hanging_pt()->master_weight(m);
                    }
                    else
                    {
                        // Loop of types of dofs
                        for (unsigned k = 0; k < n_position_type; k++)
                        {
                            // Loop over the displacement components
                            for (unsigned i = 0; i < DIM; i++)
                            {
                                position_local_eqn_at_node(k, i) = this->position_local_eqn(l, k, i);
                            }
                        }

                        // Hang weight is one
                        hang_weight = 1.0;
                    }

                    // Loop of types of dofs
                    for (unsigned k = 0; k < n_position_type; k++)
                    {
                        // Offset for faster access
                        const unsigned offset5 = dpsidxi.offset(l, k);

                        // Loop over the displacement components
                        for (unsigned i = 0; i < DIM; i++)
                        {
                            local_eqn = position_local_eqn_at_node(k, i);

                            /*IF it's not a boundary condition*/
                            if (local_eqn >= 0)
                            {
                                // Initialise the contribution
                                double sum = 0.0;

                                // Acceleration and body force
                                sum += (lambda_sq * accel[i] - b[i]) * psi(l, k);

                                // Stress term
                                for (unsigned a = 0; a < DIM; a++)
                                {
                                    unsigned count = offset5;
                                    for (unsigned b = 0; b < DIM; b++)
                                    {
                                        // Add the stress terms to the residuals
                                        sum += sigma(a, b) * interpolated_G(a, i) *
                                            dpsidxi.raw_direct_access(count);
                                        ++count;
                                    }
                                }
                                residuals[local_eqn] += W * sum * hang_weight;


                                // Get Jacobian too?
                                if (flag == 1)
                                {
                                    // Offset for faster access in general stress loop
                                    const unsigned offset1 = d_G_dX.offset(l, k, i);

                                    // Default setting for non-hanging node
                                    unsigned nn_master = 1;
                                    double hhang_weight = 1.0;

                                    // Loop over the nodes of the element again
                                    for (unsigned ll = 0; ll < n_node; ll++)
                                    {
                                        // Get pointer to local node ll
                                        Node* llocal_node_pt = this->node_pt(ll);

                                        // Cache hang status
                                        bool iis_hanging = llocal_node_pt->is_hanging();

                                        // If the node is a hanging node
                                        if (iis_hanging)
                                        {
                                            nn_master = llocal_node_pt->hanging_pt()->nmaster();
                                        }
                                        // Otherwise the node is its own master
                                        else
                                        {
                                            nn_master = 1;
                                        }


                                        // Storage for local unknown numbers at node indexed by
                                        // type and direction
                                        DenseMatrix<int> position_local_unk_at_node(n_position_type,
                                                                                     DIM);

                                        // Loop over the master nodes
                                        for (unsigned mm = 0; mm < nn_master; mm++)
                                        {
                                            if (iis_hanging)
                                            {
                                                // Find the unknown numbers
                                                position_local_unk_at_node = this->local_position_hang_eqn(
                                                    llocal_node_pt->hanging_pt()->master_node_pt(mm));

                                                // Find the hanging node weight
                                                hhang_weight =
                                                    llocal_node_pt->hanging_pt()->master_weight(mm);
                                            }
                                            else
                                            {
                                                // Loop of types of dofs
                                                for (unsigned kk = 0; kk < n_position_type; kk++)
                                                {
                                                    // Loop over the displacement components
                                                    for (unsigned ii = 0; ii < DIM; ii++)
                                                    {
                                                        position_local_unk_at_node(kk, ii) =
                                                            this->position_local_eqn(ll, kk, ii);
                                                    }
                                                }

                                                // Hang weight is one
                                                hhang_weight = 1.0;
                                            }


                                            // Loop of types of dofs again
                                            for (unsigned kk = 0; kk < n_position_type; kk++)
                                            {
                                                // Loop over the displacement components again
                                                for (unsigned ii = 0; ii < DIM; ii++)
                                                {
                                                    // Get the number of the unknown
                                                    int local_unknown =
                                                        position_local_unk_at_node(kk, ii);


                                                    /*IF it's not a boundary condition*/
                                                    if (local_unknown >= 0)
                                                    {
                                                        // Offset for faster access in general stress loop
                                                        const unsigned offset2 = d_G_dX.offset(ll, kk, ii);
                                                        const unsigned offset4 = dpsidxi.offset(ll, kk);


                                                        // General stress term
                                                        //--------------------
                                                        double sum = 0.0;
                                                        unsigned count1 = offset1;
                                                        for (unsigned a = 0; a < DIM; a++)
                                                        {
                                                            // Bump up direct access because we're only
                                                            // accessing upper triangle
                                                            count1 += a;
                                                            for (unsigned b = a; b < DIM; b++)
                                                            {
                                                                double factor =
                                                                    d_G_dX.raw_direct_access(count1);
                                                                if (a == b) factor *= 0.5;

                                                                // Offset for faster access
                                                                unsigned offset3 = d_stress_dG.offset(a, b);
                                                                unsigned count2 = offset2;
                                                                unsigned count3 = offset3;

                                                                for (unsigned aa = 0; aa < DIM; aa++)
                                                                {
                                                                    // Bump up direct access because we're only
                                                                    // accessing upper triangle
                                                                    count2 += aa;
                                                                    count3 += aa;

                                                                    // Only upper half of derivatives w.r.t.
                                                                    // symm tensor
                                                                    for (unsigned bb = aa; bb < DIM; bb++)
                                                                    {
                                                                        sum +=
                                                                            factor *
                                                                            d_stress_dG.raw_direct_access(count3) *
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
                                                        jacobian(local_eqn, local_unknown) +=
                                                            sum * W * hang_weight * hhang_weight;

                                                        // Only upper triangle (no separate test for bc as
                                                        // local_eqn is already nonnegative)
                                                        if ((i == ii) && (local_unknown >= local_eqn))
                                                        {
                                                            // Initialise the contribution
                                                            double sum = 0.0;

                                                            // Inertia term
                                                            sum += lambda_sq * time_factor * psi(ll, kk) *
                                                                psi(l, k);
                                                            sum += lambda_sq * this->dissipation_
                                                                * time_factor1 * psi(ll,kk) * psi(l,k); //TW

                                                            // Stress term
                                                            unsigned count4 = offset4;
                                                            for (unsigned a = 0; a < DIM; a++)
                                                            {
                                                                // Cache term
                                                                const double factor =
                                                                    dpsidxi.raw_direct_access(count4); // ll ,kk
                                                                ++count4;

                                                                unsigned count5 = offset5;
                                                                for (unsigned b = 0; b < DIM; b++)
                                                                {
                                                                    sum +=
                                                                        sigma(a, b) * factor *
                                                                        dpsidxi.raw_direct_access(count5); // l  ,k
                                                                    ++count5;
                                                                }
                                                            }

                                                            // Multiply by weights to form contribution
                                                            double sym_entry =
                                                                sum * W * hang_weight * hhang_weight;
                                                            // Add contribution to jacobian
                                                            jacobian(local_eqn, local_unknown) += sym_entry;
                                                            // Add to lower triangular entries
                                                            if (local_eqn != local_unknown)
                                                            {
                                                                jacobian(local_unknown, local_eqn) += sym_entry;
                                                            }
                                                        }
                                                    } // End of if not boundary condition
                                                }
                                            }
                                        }
                                    }
                                }

                            } // End of if not boundary condition

                        } // End of loop over coordinate directions
                    } // End of loop over type of dof
                } // End of loop over master nodes
            } // End of loop over nodes
        } // End of loop over integration points
    }

    double dissipation_ = 0.0; //TW

public:

    void setDissipation(double dissipation) {
        dissipation_ = dissipation;
    }
}; // end class

//==============================================================
/// FaceGeometry of the 2D RefineableQPVDElement elements
//==============================================================
template<unsigned NNODE_1D>
class FaceGeometry<RefineableQDPVDElement<2, NNODE_1D>>
    : public virtual SolidQElement<1, NNODE_1D>
{
public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<1, NNODE_1D>() {}
};

//==============================================================
/// FaceGeometry of the FaceGeometry of the 2D RefineableQPVDElement
//==============================================================
template<unsigned NNODE_1D>
class FaceGeometry<FaceGeometry<RefineableQDPVDElement<2, NNODE_1D>>>
    : public virtual PointElement
{
public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
};


//==============================================================
/// FaceGeometry of the 3D RefineableQPVDElement elements
//==============================================================
template<unsigned NNODE_1D>
class FaceGeometry<RefineableQDPVDElement<3, NNODE_1D>>
    : public virtual SolidQElement<2, NNODE_1D>
{
public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<2, NNODE_1D>() {}
};

//==============================================================
/// FaceGeometry of the FaceGeometry of the 3D RefineableQPVDElement
//==============================================================
template<unsigned NNODE_1D>
class FaceGeometry<FaceGeometry<RefineableQDPVDElement<3, NNODE_1D>>>
    : public virtual SolidQElement<1, NNODE_1D>
{
public:
    // Make sure that we call the constructor of the SolidQElement
    // Only the Intel compiler seems to need this!
    FaceGeometry() : SolidQElement<1, NNODE_1D>() {}
};

} // end namespace oomph

#endif//MERCURYDPM_REFINEABLEQDPVDELEMENT_H
