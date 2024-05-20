//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
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

#ifndef MERCURYDPM_FORMULA_HELPERS_H
#define MERCURYDPM_FORMULA_HELPERS_H

#include "GeneralDefine.h"
#include "Math/ExtendedMath.h"

namespace helpers
{
    /*!
     * \brief return type specifically for fuctions returning k and disp at once
     */
    struct KAndDisp
    {
        Mdouble k;
        Mdouble disp;
    };

    /*!
     * \brief return type specifically for fuctions returning k, disp, kt, dispt at once
     */
    struct KAndDispAndKtAndDispt
    {
        Mdouble k;
        Mdouble disp;
        Mdouble kt;
        Mdouble dispt;
    };

    /*!
     * \brief Set disp, k, dispt and kt such that is matches a given collision time tc and a normal and tangential
     * restitution coefficient r, beta for a collision of effective/reduced mass m. From Deen...Kuipers2006,
     * eq. 43 and 30
     * \todo what should be used instead of this function?
     */
    MERCURYDPM_DEPRECATED
    KAndDispAndKtAndDispt
    computeDisptFromCollisionTimeAndRestitutionCoefficientAndTangentialRestitutionCoefficientAndEffectiveMass(
            Mdouble tc, Mdouble r, Mdouble beta, Mdouble mass);

    /*!
     * \brief Calculates the maximum relative velocity allowed for a normal collision of two particles of radius r and
     * particle mass m (for higher velocities particles could pass through each other)
     * \todo what should be used instead of this function?
     */
    MERCURYDPM_DEPRECATED
    Mdouble getMaximumVelocity(Mdouble k, Mdouble disp, Mdouble radius, Mdouble mass);

    /*!
     * \brief Calculates the effective mass of a particle pair, i.e. half the harmonic mean of two particle masses.
     */
    Mdouble getEffectiveMass(Mdouble mass0, Mdouble mass1);

    /*!
     * \brief Returns the Rayleigh time step for a Hertz contact law.
     */
    Mdouble getRayleighTime(Mdouble radius, Mdouble shearModulus, Mdouble poisson, Mdouble density);

    /*!
     * \brief Returns the correct saveCount if the total number of saves, the final time and the time step is known
     */
    unsigned int getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(unsigned int numberOfSaves, Mdouble timeMax, Mdouble timeStep);

    /*!
     * \brief Returns the shear modulus calculated from Young's modulus and Poisson's ratio. 
     */
    Mdouble getShearModulus(Mdouble youngsModulus, Mdouble poisson);
}

#endif // MERCURYDPM_FORMULA_HELPERS_H
