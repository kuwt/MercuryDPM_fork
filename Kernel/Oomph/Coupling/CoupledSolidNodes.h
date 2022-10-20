//
// Created by Thomas Weinhart on 29/04/2022.
//

#ifndef MERCURY_COUPLEDSOLIDNODES_H
#define MERCURY_COUPLEDSOLIDNODES_H
#include "nodes.h"

namespace oomph
{

/**
 * adds a coupling weight and coupling force to each SolidNode
 */
class CoupledSolidNode : public SolidNode
{
private:
    /// \short Weighting factor for coupling with other method
    double coupling_weight;

    /// \short Coupling force on the node
    Vector<double> coupling_force;

public:

    CoupledSolidNode(const unsigned& n_lagrangian,
                      const unsigned& n_lagrangian_type,
                      const unsigned& n_dim,
                      const unsigned& n_position_type,
                      const unsigned& initial_n_value)
        : SolidNode(n_lagrangian, n_lagrangian_type, n_dim, n_position_type, initial_n_value) {}

    CoupledSolidNode(TimeStepper* const& time_stepper_pt_,
                      const unsigned& n_lagrangian,
                      const unsigned& n_lagrangian_type,
                      const unsigned& n_dim,
                      const unsigned& n_position_type,
                      const unsigned& initial_n_value)
        : SolidNode(time_stepper_pt_, n_lagrangian, n_lagrangian_type, n_dim, n_position_type, initial_n_value) {}

    /// \short Set and get function for the coupling weight and force
    inline void set_coupling_weight( const double& weight ) { coupling_weight = weight; }

    inline const double get_coupling_weight() { return coupling_weight; }

    inline void set_coupling_force( const Vector<double>& cForce ) { coupling_force = cForce; }

    inline Vector<double> get_coupling_force() const { return coupling_force; }
};

}// namespace oomph

#endif//MERCURY_SOLIDNODES_H
