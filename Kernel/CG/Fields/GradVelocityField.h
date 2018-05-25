//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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
#ifndef GradVelocityField_H
#define GradVelocityField_H

#include <Math/Matrix.h>
#include <Math/MatrixSymmetric.h>
#include <CG/Functions/IntegralType.h>

class BaseParticle;

class BaseInteraction;

namespace CGFields
{
    
    class GradVelocityField
    {
    public:
        
        /*!
         * \brief Default constructor, sets all field values to zero.
         */
        GradVelocityField();
        
        /*!
         * \brief Default copy constructor, copies the values of all fields.
         */
        GradVelocityField(const GradVelocityField& P) = default;
        
        /*!
         * \brief Destructor, it simply destructs the GradVelocityField and all the objects
         * it contains.
         */
        ~GradVelocityField() = default;
        
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
        GradVelocityField getSquared() const;
        
        /*!
         * \brief Copies all field values.
         */
        GradVelocityField& operator=(const GradVelocityField& P);
        
        /*!
         * \brief Adds the field values on the RHS to the LHS of the equation.
         */
        GradVelocityField& operator+=(const GradVelocityField& P);
        
        /*!
         * \brief Subtracts the field values on the RHS from the LHS of the equation.
         */
        GradVelocityField& operator-=(const GradVelocityField& P);
        
        /*!
         * \brief Divides the field values on the LHS by the RHS of the equation.
         */
        GradVelocityField& operator/=(Mdouble a);
        
        /*!
         * \brief Multiplies the field values on the left of the '*' by the
         * scalar value on the right of the '*' and returns the answer.
         */
        GradVelocityField operator*(Mdouble a) const;
        
        /*!
         * \brief This function should be called from within a loop over all
         * particles to compute all the fields that are defined as a sum over all
         * particles (e.g. density, momentum).
         */
        void addParticleStatistics(Mdouble phi, const GradVelocityField& currentInteraction);
        
        void addParticleDifferentialStatistics(Vec3D& dphi, const GradVelocityField& currentInteraction);
        
        /*!
         * \brief This function should be called from within a loop over all
         * Interactions to compute all the fields that are defined as a sum over all
         * Interactions (e.g. stress).
         */
        void addInteractionStatistics(Mdouble psi, const GradVelocityField& currentInteraction);
        
        /*!
         * \brief This function should be called from within a loop over all
         * Interactions to compute all the fields that are defined as a sum over all
         * Interactions with external objects (e.g. IFD).
         */
        void addContactPointStatistics(Mdouble phi, const GradVelocityField& currentInteraction);
        
        void setFields(const BaseInteraction& c, IntegralType type);
        
        void setCylindricalFields(const BaseInteraction& c, IntegralType type);
        
        void setFields(const BaseParticle& p);
        
        void setCylindricalFields(const BaseParticle& p);
        
        
        /*!
         * \brief Returns true if the class contains fields that are defined as a
         * sum over all Interactions (e.g. stress), else returns false.
         */
        static bool doInteractionStatistics();
        
        Mdouble getDensity() const
        {
            return density_;
        }
        
        Vec3D getMomentum() const
        {
            return momentum_;
        }
        
        Vec3D getDDensity() const
        {
            return ddensity_;
        }
        
        Matrix3D getDMomentum() const
        {
            return dmomentum_;
        }
        
        static bool evaluateFixedParticles()
        {
            return false;
        }
        
        /*!
         * A bool that determines if the derivative of the CG function has
         * to be computed
         */
        static bool isDifferentialField()
        {
            return true;
        }
    
    private:
        
        Mdouble density_;
        Vec3D momentum_;
        Vec3D ddensity_;
        Matrix3D dmomentum_;
    };
    
}
#endif
