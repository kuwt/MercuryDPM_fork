
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

#ifndef MECURYDPM_HELPERS_H
#define MECURYDPM_HELPERS_H

/*!
 * \todo gmb Remove this.
 *
 * This inclusion is here to ensure backwards compatibility.
 * Some drivers needs to be modified to remove this todo.
 */
#include "Helpers/Helpers.h"

/*!
 * \todo gmb Remove these? If not move them to Helpers/FormulaHelpers.h
 */
//    /*!
//     * \brief Set disp and k such that is matches a given collision time tc and restitution coefficient r
//     * for a collision of effective/reduced mass m.
//     * \deprecated use species->setCollisionTimeAndRestitutionCoefficient
//     *   (collisionTime, dissipationTimeScale, 2.0*effectiveMass) instead
//     */
//    MERCURYDPM_DEPRECATED
//    KAndDisp computeKAndDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the collision time for a given stiffness, dissipation, and effective mass
//     * \deprecated use species->computeCollisionTime(2.0*effectiveMass) instead
//     * \todo This does not result in the same value as the given alternative.
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeCollisionTimeFromKAndDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass);
//
//    /*!
//     * \brief Calculates the restitution coefficient time for a given stiffness, dissipation, and effective mass
//     * \deprecated use species->computeRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeRestitutionCoefficientFromKAndDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass);
//
//    /*!
//     * \brief Calculates the dissipation for a given stiffness, restitution coefficient, and effective mass
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeDispFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the collision time for a given stiffness, restitution coefficient, and effective mass
//     * \deprecated use species->computeCollisionTime(2.0*effectiveMass) instead
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeCollisionTimeFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the dissipation for a given stiffness, collision time, and effective mass
//     * \deprecated use species->setStiffnessAndRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeDispFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass);
//
//    /*!
//     * \brief Calculates the restitution coefficient for a given stiffness, collision time, and effective mass
//     * \deprecated use species->computeRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeRestitutionCoefficientFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass);
//
//    /*!
//     * \brief Calculates the dissipation for a given collision time, restitution coefficient, and effective mass
//     * \deprecated use species->setCollisionTimeAndRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the stiffness for a given collision time, restitution coefficient, and effective mass
//     * \deprecated use species->setCollisionTimeAndRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeKFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the stiffness for a given collision time, dissipation, and effective mass
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeKFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass);
//
//    /*!
//     * \brief Calculates the resitution coefficient for a given collision time, dissipation, and effective mass
//     * \deprecated use species->computeRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeRestitutionCoefficientFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass);
//
//    /*!
//     * \brief Calculates the stiffness for a given dissipation, restitution coefficient, and effective mass
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeKFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the collision time for a given dissipation, restitution coefficient, and effective mass
//     * \deprecated use species->computeCollisionTime(2.0*effectiveMass) instead
//     */
//    MERCURYDPM_DEPRECATED
//    Mdouble computeCollisionTimeFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass);

#endif