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

#ifndef SuperQuad_H
#define SuperQuad_H

#include "Math/SmallVector.h"
#include "BaseParticle.h"

typedef Vec3D BodyFixedCoordinates;
typedef Vec3D LabFixedCoordinates;

/*!
 * \class SuperQuad
 * \brief Class that implements superquadric particles, which are non-spherical
 */
class SuperQuadric : public BaseParticle
{
public:
    /*!
     * \brief Basic Particle constructor, creates an Particle at (0,0,0) with radius, mass and inertia equal to 1
     */
    SuperQuadric();

    /*!
     * \brief Particle copy constructor, which accepts as input a reference to a Particle. It creates a copy of this Particle and all it's information. Usually it is better to use the copy() function for polymorfism.
     */
    SuperQuadric(const SuperQuadric& p);

    /*!
     * \brief Particle destructor, needs to be implemented and checked if it removes tangential spring information
     */
    ~SuperQuadric() override;

    /*!
     * \brief Particle copy method. It calls to copy constructor of this Particle, useful for polymorfism
     */
    SuperQuadric* copy() const override;
    
    void write(std::ostream& os) const override;
	
    void read(std::istream& is) override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;

    /*!
     * \brief We consider the superellipsoid definition stated in Chapter 2 of the book segmentation and recovery of superquadrics by Jaklic and co.
     */

    void setAxesAndExponents(const Mdouble& a1, const Mdouble& a2, const Mdouble& a3, const Mdouble& eps1, const Mdouble& eps2);
    
    void setAxesAndExponents(const Vec3D& axes, const Mdouble& eps1, const Mdouble& eps2);

    void setAxes(const Mdouble& a1, const Mdouble& a2, const Mdouble& a3);

    void setAxes(const Vec3D& axes);

    void setExponents(const Mdouble& eps1, const Mdouble& eps2);

    ///\todo TW we could remove this function from the BaseParticle and use a dynamic_cast instead
    Vec3D getAxes() override;

    Mdouble getExponentEps1() override;

    Mdouble getExponentEps2() override;
    
    Mdouble getVolume() const override;
    
    void setInertia() override;
    
    ///\brief return the radius of the sphere that fits precisely around the particle. Currently only implemented for
    ///ellipsoids
    Mdouble getInteractionRadius() const override;
    
    std::vector<BaseInteraction *> getInteractionWith(BaseParticle *const P,
                                                                                         const unsigned timeStamp,
                                                                                         InteractionHandler *const interactionHandler) override;
    
    Mdouble getCurvature(const LabFixedCoordinates& xGlobal) const;
    
    bool isInContactWith(const BaseParticle* const p) const override;
    
    SmallVector<3> computeShapeGradientLocal(const LabFixedCoordinates& xGlobal) const;
    
private:
    
    ///returns the interaction between two superquads: if they do not touch, it returns nullptr, otherwise it returns the
    ///pointer to the interaction between this particle and the given pQuad.
    std::vector<BaseInteraction*> getInteractionWithSuperQuad(SuperQuadric* const pQuad, const unsigned timeStamp,
                                                 InteractionHandler* const interactionHandler);
    
    
    SmallMatrix<3,3> computeHessian(const LabFixedCoordinates& xGlobal) const;
    
    Mdouble computeShape(const LabFixedCoordinates& xGlobal) const;
    
    SmallVector<4> functionThatShouldBecomeZeroForContactDetection(const SmallVector<4>& approximateContact,
                                                                       const SuperQuadric* const pQuad) const;
    
    SmallMatrix<4, 4> computeJacobian(const SmallVector<4>& approximateContact, const SuperQuadric* const pQuad) const;
    
    SmallVector<4> getInitialGuessForContact(const SuperQuadric* pQuad, BaseInteraction* const C) const;
    
    Mdouble computeOverlapAlpha(const LabFixedCoordinates& contactPoint, const LabFixedCoordinates& normal) const;
    
    SmallVector<4> getContactPoint(const SuperQuadric* const pQuad, BaseInteraction* C) const;
    
    Mdouble eps1_,eps2_;
    
    Vec3D axes_; /// the three major or minor axes (a1,a2,a3)
};
#endif
