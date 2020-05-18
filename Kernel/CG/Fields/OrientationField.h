//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef MERCURY_ORIENTATIONFIELD_H
#define MERCURY_ORIENTATIONFIELD_H

#include <Math/Matrix.h>
#include <Math/MatrixSymmetric.h>
#include <CG/Functions/IntegralType.h>
#include <array>

class BaseParticle;

class BaseInteraction;
namespace CGFields
{

/*!
 * \brief Contains the computed field values, like density, momentum and stress.
 * \details CGPoints inherits from this class; CGPoints::evaluate adds to the
 * values of these variables.
 * \todo These are currently the only fields that are computed.
 * However, this class is destined to be extended to
 * contain additional information such as fabric, energy, local angular
 * momentum. Also, a simpler version is planned, where only particle statistics
 * are evaluated (density and momentum).
 */
    
    class OrientationField
    {
    public:
        OrientationField();
        
        OrientationField(const OrientationField& other) = default;
        
        ~OrientationField() = default;
        
        static void writeNames(std::ostream& os, unsigned countVariables);
        
        /*!
         * \brief Writes class content into an output stream, typically a stat file.
         */
        void write(std::ostream& os) const;
        
        /*!
         * \brief Writes human-readable class content into an output stream, typically a stat file.
         */
        void output(std::ostream& os) const;
        
        /*!
         * \brief Sets all fields to zero.
         */
        void setZero();
        
        /*!
         * \brief Returns the square of all field values (to calculate standard deviation).
        */
        OrientationField getSquared() const;
        
        /*!
         * \brief Copies all field values.
         */
        OrientationField& operator=(const OrientationField& P);
        
        /*!
         * \brief Adds the field values on the RHS to the LHS of the equation.
         */
        OrientationField& operator+=(const OrientationField& P);
        
        /*!
         * \brief Subtracts the field values on the RHS from the LHS of the equation.
         */
        OrientationField& operator-=(const OrientationField& P);
        
        /*!
         * \brief Divides the field values on the LHS by the RHS of the equation.
         */
        OrientationField& operator/=(Mdouble a);
        
        /*!
         * \brief Multiplies the field values on the left of the '*' by the
         * scalar value on the right of the '*' and returns the answer.
         */
        OrientationField operator*(Mdouble a) const;
        
        /*!
         * \brief This function should be called from within a loop over all
         * particles to compute all the fields that are defined as a sum over all
         * particles (e.g. density, momentum).
         */
        void addParticleStatistics(Mdouble phi, const OrientationField& currentInteraction);
        
        void setFields(const BaseParticle& p);
        
        void setCylindricalFields(const BaseParticle& p);
        
        
        MatrixSymmetric3D getOrientation() const
        {
            return orientation_;
        }
        
        
        static bool evaluateFixedParticles()
        {
            return false;
        }
        
        static bool doInteractionStatistics()
        {
            return false;
        }
        
        /*!
         * A bool that determines if the derivative of the CG function has
         * to be computed
         */
        static bool isDifferentialField()
        {
            return false;
        }
        
        void setFields(const BaseInteraction& c, IntegralType type)
        {}
        
        void setCylindricalFields(const BaseInteraction& c, IntegralType type)
        {}
        
        void addParticleDifferentialStatistics(Vec3D& dphi, const OrientationField& currentInteraction)
        {}
        
        void addInteractionStatistics(Mdouble psi, const OrientationField& currentInteraction)
        {}
        
        void addContactPointStatistics(Mdouble phi, const OrientationField& currentInteraction)
        {}
    
    private:
        MatrixSymmetric3D orientation_;
    };
}

#endif //MERCURY_ORIENTATIONFIELD_H
