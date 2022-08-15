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

#ifndef MERCURY_ANISOTROPICHOOKEAN_H
#define MERCURY_ANISOTROPICHOOKEAN_H
#include "constitutive_laws.h"

namespace oomph
{

/**
 * Adds anisotropy to the general Hookean constitutive law
 */
class AnisotropicHookean : public oomph::GeneralisedHookean
{
public:
    /// The constructor takes the pointers to values of material parameters:
    /// Poisson's ratio and Young's modulus.
    AnisotropicHookean( double* nu_pt, double* e_pt )
        : GeneralisedHookean( nu_pt, e_pt ) {}

    //=====================================================================
    /// Calculate the contravariant 2nd Piola Kirchhoff
    /// stress tensor. Arguments are the
    /// covariant undeformed (stress-free) and deformed metric
    /// tensors, g and G, and the matrix in which to return the stress tensor.
    //=====================================================================
    void calculate_second_piola_kirchhoff_stress(
        const DenseMatrix<double>& g,
        const DenseMatrix<double>& G,
        DenseMatrix<double>& sigma)
    {
        GeneralisedHookean::calculate_second_piola_kirchhoff_stress(g, G, sigma);

        // Make anisotropic
        const unsigned dim = sigma.nrow();
        for (unsigned i = 0; i < dim; i++)
        {
            for ( unsigned j = 0; j < dim; j++ )
            {
                sigma( i, j) *= anisotropy[i];
            }
        }

        // Symmetrize
        for (unsigned i = 0; i < dim; i++)
        {
            for (unsigned j = 0; j < i; j++)
            {
                sigma(i, j) = sigma(j, i);
            }
        }
    }

    std::array<double, 3> anisotropy { 1.0, 1.0, 1.0 };
};
}

#endif//MERCURY_ANISOTROPICHOOKEAN_H
