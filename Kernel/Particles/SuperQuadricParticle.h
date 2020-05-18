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

#ifndef SuperQuad_H
#define SuperQuad_H

#include "Math/SmallVector.h"
#include "BaseParticle.h"
#include "BaseInteractable.h"
#include "Species/ParticleSpecies.h"

typedef Vec3D BodyFixedCoordinates;
typedef Vec3D LabFixedCoordinates;

/*!
 * \class SuperQuad
 * \brief Class that implements superquadric particles, which are non-spherical.
 *
 * \details We use the super-ellipsoid definition stated in Chapter 2 of the book "Segmentation and recovery of
 * superquadrics" by Jaklic et al. This class contains the geometrical information of the particles, i.e. the length of
 * the principal axes and the "blockiness" parameters.
 *
 * This class computes the volume and inertia tensor of the particle, and also the contact-point between two
 * superquadrics, or a normal particle and superquadric. It also computes its interaction-radius and curvature at a
 * certain position. It should be noted that two coordinate systems are used in the methods of this class: the
 * lab-fixed coordinates, which are points in physical space as seen by an outside observer, and the body-fixed
 * coordinates, which are points in space relative to the (not-oriented) superquadric particle.
 *
 * Most of the contact-detection algorithm is based on Comp. Part. Mech. (2017) 4 : 101-118, "Efficient implementation
 * of superquadric particles in Discrete element Method within an open-source framework" by Podlozhnyuk, Pirker and
 * Kloss.
 */
class SuperQuadricParticle final : public BaseParticle
{
public:
    /*!
     * \brief Basic Particle constructor, creates a superquadric with axes (1,1,1) and exponents (2,2), so it creates a
     * sphere with radius 1.
     */
    SuperQuadricParticle();
    
    /*!
     * \brief Copy constructor, which accepts as input a reference to a Superquadric.
     * It creates a copy of this Particle and all it's information.
     * Usually it is better to use the copy() function for polymorphism.
     */
    SuperQuadricParticle(const SuperQuadricParticle& p);
    
    SuperQuadricParticle(const BaseParticle& p);
    
    /*!
     * \brief Destructor, needs to be implemented and checked to see if it is the largest or smallest particle currently
     * in its particleHandler
     */
    ~SuperQuadricParticle() override;
    
    /*!
     * \brief Copy method. It calls to copy constructor of this superquadric, useful for polymorphism
     */
    SuperQuadricParticle* copy() const override;
    
    /*!
     * \brief Write function: write this superquadric to the given output-stream, for example a restart-file
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Read function: read in the information for this superquadric from the given input-stream, for example a
     * restart file
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Returns the name of the class, here "SuperQuadricParticle"
     */
    std::string getName() const override;
    
    /*!
     * \brief Set the geometrical properties of the superquadrics, namely the axes-lengths a1, a2 and a3, and the
     * exponents epsilon1 and epsilon2.
     * We use the super-ellipsoid definition stated in Chapter 2 of the book "Segmentation and recovery of
     * superquadrics" by Jaklic et al.
     */
    void setAxesAndExponents(const Mdouble& a1, const Mdouble& a2, const Mdouble& a3, const Mdouble& eps1,
                             const Mdouble& eps2);
    
    /*!
     * \brief Set the geometrical properties of the superquadrics, namely the axes-lengths axes, and the
     * exponents epsilon1 and epsilon2.
     * We use the super-ellipsoid definition stated in Chapter 2 of the book "Segmentation and recovery of
     * superquadrics" by Jaklic et al.
     */
    void setAxesAndExponents(const Vec3D& axes, const Mdouble& eps1, const Mdouble& eps2);
    
    /*!
     * \brief Set the axes-lengths to a1, a2 and a3 for this superquadric.
     * We use the super-ellipsoid definition stated in Chapter 2 of the book "Segmentation and recovery of
     * superquadrics" by Jaklic et al.
     */
    void setAxes(const Mdouble& a1, const Mdouble& a2, const Mdouble& a3);
    
    /*!
     * \brief Set the axes-lengths to axes for this superquadric.
     * We use the super-ellipsoid definition stated in Chapter 2 of the book "Segmentation and recovery of
     * superquadrics" by Jaklic et al.
     */
    void setAxes(const Vec3D& axes) override;
    
    /*!
     * \brief Set the exponents to eps1 and eps2 for this superquadric.
     * We use the super-ellipsoid definition stated in Chapter 2 of the book "Segmentation and recovery of
     * superquadrics" by Jaklic et al.
     */
    void setExponents(const Mdouble& eps1, const Mdouble& eps2) override;
    
    ///\todo TW we could remove this function from the BaseParticle and use a dynamic_cast instead
    ///\todo ID Middle-term plan is to template the BaseParticle on shape-type, so that we won't have to cast etc.
    /*!
     * \brief Get the axes-lengths of this superquadric.
     * We use the super-ellipsoid definition stated in Chapter 2 of the book "Segmentation and recovery of
     * superquadrics" by Jaklic et al.
     */
    Vec3D getAxes() const override;
    
    /*!
 * \brief Get the first exponent of this superquadric.
     * We use the super-ellipsoid definition stated in Chapter 2 of the book "Segmentation and recovery of
     * superquadrics" by Jaklic et al.
 */
    Mdouble getExponentEps1() const override;
    
    /*!
 * \brief Get the second exponent of this superquadric.
     * We use the super-ellipsoid definition stated in Chapter 2 of the book "Segmentation and recovery of
     * superquadrics" by Jaklic et al.
 */
    Mdouble getExponentEps2() const override;
    
    /*!
     * Get the volume of this superquadric.
     */
    Mdouble getVolume() const override;
    
    /*!
     * \brief Compute and set the inertia-tensor for this superquadric. For internal use only.
     */
    void setInertia() override;
    
    void setRadius(const Mdouble radius) override;
    
    /*!
     * \brief Checks if this superquadric is in interaction with the given particle, and if
     * so, returns vector of pointer to the associated BaseInteraction object (else returns empty vector).
     */
    BaseInteraction* getInteractionWith(BaseParticle* P, unsigned timeStamp,
                                                     InteractionHandler* interactionHandler) override;
    
    /*!
     * \brief Get the mean curvature of this superquadric at the given (lab-fixed) position, see Podlozhyuk et al.
     * (2017) eq (39)
     */
    Mdouble getCurvature(const LabFixedCoordinates& labFixedCoordinates) const override;
    
    /*!
     * \brief Get whether or not this superquadric is in contact with the given particle.
     */
    bool isInContactWith(const BaseParticle* p) const override;
    
    /*!
     * \brief Compute and get the gradient of the shape-function at the given (lab-fixed) position.
     */
    SmallVector<3> computeShapeGradientLabFixed(const LabFixedCoordinates& labFixedCoordinates) const;
    
    
    /*!
     * \brief Checks if this superquadric is in interaction with the given superquadric, and if
     * so, returns vector of pointer to the associated BaseInteraction object (else returns empty vector).
     */
    BaseInteraction* getInteractionWithSuperQuad(SuperQuadricParticle* p, unsigned timeStamp,
                                                              InteractionHandler* interactionHandler);
    
    /*!
     * \brief Compute and get the hessian ("second derivative") of the shape-function at the given (lab-fixed) position.
     */
    SmallMatrix<3, 3> computeHessianLabFixed(const LabFixedCoordinates& labFixedCoordinates) const;
    
    /*!
     * \brief Compute and get the shape-functiion at the given (lab-fixed) position.
     */
    Mdouble computeShape(const LabFixedCoordinates& labFixedCoordinates) const;
    
    /*!
     * \brief Objective function for contact detection between the two given superquadrics. See  Podlozhyuk et al.
     * (2017) eq (22).
     */
    SmallVector<4> computeResidualContactDetection(const SmallVector<4>& position,
                                                   const SuperQuadricParticle* p1,
                                                   const SuperQuadricParticle* p2) const;
    
    /*!
     * \brief Compute and return the derivative of functionThatShouldBecomeZeroForContactDetection, both to the position
     * and the Lagrange multiplier, and evaluated at the contact point.
     */
    SmallMatrix<4, 4> getJacobianOfContactDetectionObjective(const SmallVector<4>& contactPoint,
                                                             const SuperQuadricParticle* p1,
                                                             const SuperQuadricParticle* p2) const;
    
    /*!
     * \brief Get an initial guess for the contact-point between this particle and the given particle.
     */
    SmallVector<4> getInitialGuessForContact(const SuperQuadricParticle* pQuad, BaseInteraction* C) const;
    
    /*!
     * \brief Compute the distance between the contact-point and surface of this superquadric particle.
     */
    Mdouble overlapFromContactPoint(const LabFixedCoordinates& contactPoint, const LabFixedCoordinates& normal) const;
    
    /*!
     * \brief Compute the contact point between this and the given superquadric particle.
     */
    SmallVector<4> getContactPoint(const SuperQuadricParticle* p, BaseInteraction* C) const;
    
    /*!
     * \brief If the "normal" procedure fails to find a contact point, use an alternative approach that involves
     * starting with two spheres to compute the interaction, and becoming less and less spherical.
     */
    SmallVector<4> getContactPointPlanB(const SuperQuadricParticle* pOther, unsigned numberOfSteps) const;
    
    /*!
     * \brief Perform the actual Newton-iterations to find the contact point. Note, that it is given back as a
     * parameter.
     */
    bool computeContactPoint(SmallVector<4>& contactPoint, const SuperQuadricParticle* p1,
                             const SuperQuadricParticle* p2) const;
    
    
    void writeDebugMessageStep1(const SuperQuadricParticle* pQuad, const SmallVector<4>& contactPointPlanB) const;
    
    void writeDebugMessageStep2(const SuperQuadricParticle* pQuad, const Vec3D& dAxesThis, const Mdouble& dn11,
                                const Mdouble& dn12, const Vec3D& dAxesOther, const Mdouble& dn21,
                                const Mdouble& dn22) const;
    
    void writeDebugMessageStep3(const Vec3D& axesThis, const Mdouble& n11, const Mdouble& n12, const Vec3D& axesOther,
                                const Mdouble& n21, const Mdouble& n22) const;
    
    void
    writeDebugMessageMiddleOfLoop(const SuperQuadricParticle& p1, const SuperQuadricParticle& p2, SmallVector<4>& contactPointPlanB,
                                  const unsigned int& counter) const;

    /**
     * \brief returns the radius plus half the interactionDistance of the mixed species
     */
    Mdouble getInteractionRadius(const BaseParticle* particle) const;

    /// Computes the particle's (inverse) mass and inertia.
    void computeMass(const ParticleSpecies& s) override;

private:
    
    /*!\brief Get the radius of the sphere that fits precisely around the particle.
     * \todo Currently only implemented for ellipsoids
     */
    void setBoundingRadius();
    
    /*!
     * \brief Blockiness parameters
     * \details Blockiness parameters should be in the range (0,1], where a sphere or ellipsoid is represented by
     * eps1_ = eps2_ = 1.
     */
    Mdouble eps1_, eps2_;
    
    /*!
     * \brief Lengths of principal axes (a1, a2, a3).
     */
    Vec3D axes_;
};

#endif
