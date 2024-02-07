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

#ifndef MERCURYDPM_DisplacementField_H
#define MERCURYDPM_DisplacementField_H

#include <Math/Matrix.h>
#include <Math/MatrixSymmetric.h>
#include <CG/Functions/IntegralType.h>
#include <array>

#include "BaseFields.h"
#include "CG/CGHandler.h"
#include "BaseHandler.h"
#include "DPMBase.h"

class DPMBase;

class BaseParticle;

class BaseInteraction;
namespace CGFields
{

/*!
 * \brief Computed the displacement fields
 * \details The displacement moment and the displacement moment flux
 * are calculated from the previous position of the particle and the
 * time since the last evaluation of the data files
 */
    
    class DisplacementField : public BaseFields
    {
    public:
        DisplacementField();
        
        DisplacementField(const DisplacementField& other) = default;
        
        ~DisplacementField() = default;
        
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
        DisplacementField getSquared() const;
        
        /*!
         * \brief Copies all field values.
         */
        DisplacementField& operator=(const DisplacementField& P);
        
        /*!
         * \brief Adds the field values on the RHS to the LHS of the equation.
         */
        DisplacementField& operator+=(const DisplacementField& P);
        
        /*!
         * \brief Subtracts the field values on the RHS from the LHS of the equation.
         */
        DisplacementField& operator-=(const DisplacementField& P);
        
        /*!
         * \brief Divides the field values on the LHS by the RHS of the equation.
         */
        DisplacementField& operator/=(Mdouble a);
        
        /*!
         * \brief Multiplies the field values on the left of the '*' by the
         * scalar value on the right of the '*' and returns the answer.
         */
        DisplacementField operator*(Mdouble a) const;
        
        /*!
         * \brief This function should be called from within a loop over all
         * particles to compute all the fields that are defined as a sum over all
         * particles (e.g. density, momentum).
         */
        void addParticleStatistics(Mdouble phi, const DisplacementField& currentInteraction);
        
        void setFields(const BaseParticle& p);
        
        void setCylindricalFields(const BaseParticle& p);
        
        
        MatrixSymmetric3D getDisplacementMomentumFlux() const { return displacementMomentumFlux_; }
        Vec3D getDisplacementMomentum() const { return displacementMomentum_; }
        
        
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
        
        void addParticleDifferentialStatistics(Vec3D& dphi, const DisplacementField& currentInteraction)
        {}
        
        void addInteractionStatistics(Mdouble psi, const DisplacementField& currentInteraction)
        {}
        
        void addContactPointStatistics(Mdouble phi, const DisplacementField& currentInteraction)
        {}
    
    private:
        MatrixSymmetric3D displacementMomentumFlux_;
        Vec3D displacementMomentum_;
    };
}

#endif //MERCURYDPM_DisplacementField_H
