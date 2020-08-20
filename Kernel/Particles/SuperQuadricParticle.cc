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

#include "DPMBase.h"
#include "BaseParticle.h"
#include <cmath>
#include "SuperQuadricParticle.h"
#include "InteractionHandler.h"
#include "ParticleHandler.h"
#include "SpeciesHandler.h"

/*!
 * \details calls the default constructor of BaseParticle, and creates an SuperEllipsoid with axes (1,1,1) and exponents
 * (2,2), so it creates a sphere with radius 1.
 */
SuperQuadricParticle::SuperQuadricParticle()
        : BaseParticle()
{
    axes_ = Vec3D(1.0, 1.0, 1.0);
    eps1_ = 1.0;
    eps2_ = 1.0;
    logger(DEBUG, "SuperQuadricParticle::SuperQuadricParticle() finished");
}

/*!
 * \details Constructor that copies most of the properties of the given particle.
 *          Please note that not everything is copied, for example the position 
 *          in the HGrid is not determined yet by the end of this constructor. 
 *          It also does not copy the interactions and the pointer to the handler
 *          that handles this particle. Use with care.
 * \param[in,out] p  Reference to the SuperQuad this one should become a copy of.
 */
SuperQuadricParticle::SuperQuadricParticle(const SuperQuadricParticle& p)
        : BaseParticle(p)
{
    axes_ = p.axes_;
    eps1_ = p.eps1_;
    eps2_ = p.eps2_;
}


SuperQuadricParticle::SuperQuadricParticle(const BaseParticle& p) : BaseParticle(p)
{
    Mdouble radius = p.getRadius();
    axes_ = Vec3D(radius, radius, radius);
    eps1_ = 1.0;
    eps2_ = 1.0;
}

/*!
 * \details Destructor. It asks the ParticleHandler to check if this was the
 *          smallest or largest particle and adjust itself accordingly.
 */
SuperQuadricParticle::~SuperQuadricParticle()
{
    if (getHandler() != nullptr)
    {
        getHandler()->checkExtremaOnDelete(this);
    }
    logger(DEBUG, "SuperQuadricParticle::SuperQuadricParticleParticle() of particle % finished.", getId());
    
}

/*!
 * \details Copy method. Uses copy constructor to create a copy on the heap. Useful for polymorphism.
 * \return pointer to the particle's copy
 */
SuperQuadricParticle* SuperQuadricParticle::copy() const
{
    return new SuperQuadricParticle(*this);
}

/*!
 * \details SuperQuadricParticle print method, which accepts an std::ostream as input. It prints human readable SuperQuadricParticle
 * information to the given output-stream.
 * \param[in,out] os    stream to which the info is written, e.g. a restart-file or std::cout.
 */
void SuperQuadricParticle::write(std::ostream& os) const
{
    BaseParticle::write(os);
    os << " axes " << axes_
       << " exp1 " << eps1_
       << " exp2 " << eps2_;
}

/*!
 * \details Returns the name of the object; in this case "SuperQuadricParticle".
 * \return The std::string "SuperQuadricParticle".
 */
std::string SuperQuadricParticle::getName() const
{
    return "SuperQuadricParticle";
}


/*!
 * \details Particle read function, which reads the axes_ and both epsilons from the given input-stream.
 * \param[in,out] is    input stream with particle properties, e.g. a restart-file.
 */
void SuperQuadricParticle::read(std::istream& is)
{
    BaseParticle::read(is);
    std::string dummy;
    is >> dummy >> axes_ >> dummy >> eps1_ >> dummy >> eps2_;
    logger.assert_always(eps1_ < 1 + 1e-10, "epsilon1 should be at most 1");
    logger.assert_always(eps2_ < 1 + 1e-10, "epsilon2 should be at most 1");
}

// Set and Get Functions

void SuperQuadricParticle::setAxesAndExponents(const Mdouble& a1, const Mdouble& a2, const Mdouble& a3, const Mdouble& eps1,
                                       const Mdouble& eps2)
{
    eps1_ = eps1;
    eps2_ = eps2;
    logger.assert_always(eps1_ < 1 + 1e-10, "epsilon1 should be at most 1");
    logger.assert_always(eps2_ < 1 + 1e-10, "epsilon2 should be at most 1");

    setAxes(a1,a2,a3);
}

void SuperQuadricParticle::setAxesAndExponents(const Vec3D& axes, const Mdouble& eps1, const Mdouble& eps2)
{
    eps1_ = eps1;
    eps2_ = eps2;
    logger.assert_always(eps1_ < 1 + 1e-10, "epsilon1 should be at most 1");
    logger.assert_always(eps2_ < 1 + 1e-10, "epsilon2 should be at most 1");

    setAxes(axes);
}

void SuperQuadricParticle::setAxes(const Vec3D& axes)
{
    axes_ = axes;
    setBoundingRadius();
    if (getSpecies() != nullptr)
    {
        getSpecies()->computeMass(this);
    }
}

void SuperQuadricParticle::setAxes(const Mdouble& a1, const Mdouble& a2, const Mdouble& a3)
{
    setAxes({a1,a2,a3});
}

void SuperQuadricParticle::setExponents(const Mdouble& eps1, const Mdouble& eps2)
{
    eps1_ = eps1;
    eps2_ = eps2;
    setBoundingRadius();
    logger.assert_always(eps1_ < 1 + 1e-10, "epsilon1 should be at most 1");
    logger.assert_always(eps2_ < 1 + 1e-10, "epsilon2 should be at most 1");
    if (getSpecies() != nullptr)
    {
        getSpecies()->computeMass(this);
    }
}

Vec3D SuperQuadricParticle::getAxes() const
{
    return axes_;
}

Mdouble SuperQuadricParticle::getExponentEps1() const
{
    return eps1_;
}

Mdouble SuperQuadricParticle::getExponentEps2() const
{
    return eps2_;
}

/*!
 * \details Returns the volume of the SuperEllipsoid, which is calculated using its principal axis and
 * the exponents. The analytical expressions are taken from Chapter 2 of the book "Segmentation and Recovery of
 * Superquadrics" by <a href="https://homes.di.unimi.it/borghese/Teaching/DigitalAnimation/Old/DigitalAnimation_2002_2003/00_Superquadriche.pdf"> Jaklic et al </a>.
 * However, the beta functions that are part of these expressions are approximations, see ExtendedMath.cc
 * \return The actual volume of this superquadric.
 */
Mdouble SuperQuadricParticle::getVolume() const
{
    logger.assert(getHandler() != nullptr, "[SuperQuadricParticle::getVolume] no particle handler specified");

    return (2.0 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) * mathsFunc::beta(0.5 * eps1_ + 1.0, eps1_) *
           mathsFunc::beta(0.5 * eps2_, 0.5 * eps2_);
}

/*!
 * \details This function computes the principal moments of inertia for the superellipsoids. Again, the analytical
 * expressions are taken from Chapter 2 (pg. 36) of the book "Segmentation and Recovery of Superquadrics" by by <a href="https://homes.di.unimi.it/borghese/Teaching/DigitalAnimation/Old/DigitalAnimation_2002_2003/00_Superquadriche.pdf"> Jaklic et al </a>.
 */
void SuperQuadricParticle::setInertia()
{
    MatrixSymmetric3D inertia;
    
    inertia.XX = getSpecies()->getDensity() * (0.5 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) *
                 (axes_.Y * axes_.Y * mathsFunc::beta(1.5 * eps2_, 0.5 * eps2_) *
                  mathsFunc::beta(0.5 * eps1_, 2.0 * eps1_ + 1.0)
                  + 4.0 * axes_.Z * axes_.Z * mathsFunc::beta(0.5 * eps2_, 0.5 * eps2_ + 1) *
                    mathsFunc::beta(1.5 * eps1_, eps1_ + 1.0));
    
    inertia.YY = getSpecies()->getDensity() * (0.5 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) *
                 (axes_.X * axes_.X * mathsFunc::beta(1.5 * eps2_, 0.5 * eps2_) *
                  mathsFunc::beta(0.5 * eps1_, 2.0 * eps1_ + 1.0)
                  + 4.0 * axes_.Z * axes_.Z * mathsFunc::beta(0.5 * eps2_, 0.5 * eps2_ + 1) *
                    mathsFunc::beta(1.5 * eps1_, eps1_ + 1.0));
    
    inertia.ZZ = getSpecies()->getDensity() * (0.5 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) *
                 (axes_.X * axes_.X + axes_.Y * axes_.Y) * mathsFunc::beta(1.5 * eps2_, 0.5 * eps2_) *
                 mathsFunc::beta(0.5 * eps1_, 2.0 * eps1_ + 1.0);
    
    BaseParticle::setInertia(inertia);
}

void SuperQuadricParticle::setRadius(const Mdouble radius)
{
    logger(ERROR,"This function should not be used");
}

void SuperQuadricParticle::setBoundingRadius()
{

    if(mathsFunc::isEqual(eps2_,1,std::numeric_limits<Mdouble>::epsilon()) && mathsFunc::isEqual(eps1_,1,std::numeric_limits<Mdouble>::epsilon()))
    {
        BaseParticle::setRadius(std::max(std::max(axes_.Y,axes_.X),axes_.Z));
        return;
    }else
    {

        const Mdouble axesX = std::max(axes_.X, axes_.Y);
        const Mdouble axesY = std::min(axes_.X, axes_.Y);

        Mdouble alpha;

        const Mdouble eps1 = std::min(.96, eps1_);
        const Mdouble eps2 = std::min(.96, eps2_);

        alpha = std::pow(axesY / axesX, 2.0 / (2.0 / eps2 - 2.0));

        const Mdouble help1 = std::pow(alpha, 2.0 / eps2);
        const Mdouble gamma = std::pow(1.0 + help1, eps2 / eps1 - 1.0);
        const Mdouble beta = std::pow(gamma * axes_.Z * axes_.Z / (axesX * axesX), 1.0 / (2.0 / eps1 - 2.0));
        const Mdouble xTilde = std::pow(std::pow(1 + help1, eps2 / eps1) + std::pow(beta, 2.0 / eps1),
                                        -eps1 / 2.0);
        BaseParticle::setRadius(std::sqrt(mathsFunc::square(axesX * xTilde)
                                          + mathsFunc::square(alpha * axesY * xTilde)
                                          + mathsFunc::square(beta * axes_.Z * xTilde)));

    }
}

/*!
 * \details Overwrites BaseInteractable::getInteractionWith. First checks if the bounding radii overlap, and if so,
 * changes the given particle into a superquadric and calls getInteractionWithSuperQuadric to obtain the relevant
 * contact data.
 * \param p                     Pointer to the particle we want to get the interaction with, can be any type of particle
 * \param timeStamp             Time stamp to be assigned to the interaction object (i.e., the current time)
 * \param interactionHandler    BaseInteraction container from where the interaction is retrieved, and to which it is
 *                              assigned (if it is a new interaction).
 * \return                      A vector with size 1 (if there is an interaction) or 0 (if there is no interaction).
 */
BaseInteraction* SuperQuadricParticle::getInteractionWith(BaseParticle* const p, const unsigned timeStamp,
                                                               InteractionHandler* const interactionHandler)
{
    //get the normal (from P away from the contact)
    const LabFixedCoordinates branchVector = p->getPosition() - getPosition();
    //Get the square of the distance between particle i and particle j
    const Mdouble distanceSquared = Vec3D::getLengthSquared(branchVector);
    auto mixedSpecies = getSpecies()->getHandler()->getMixedObject(getSpecies(),p->getSpecies());

    const Mdouble sumOfInteractionRadii = getMaxInteractionRadius()+p->getMaxInteractionRadius();
    if (distanceSquared < (sumOfInteractionRadii * sumOfInteractionRadii))
    {
        if (p->isSphericalParticle())
        {
            SuperQuadricParticle* pQuad = new SuperQuadricParticle;
            pQuad->setAxes(p->getRadius(), p->getRadius(), p->getRadius());
            BaseInteraction* contacts = getInteractionWithSuperQuad(pQuad, timeStamp, interactionHandler);
            delete pQuad;
            return contacts;
        } else {
            SuperQuadricParticle* pQuad = static_cast<SuperQuadricParticle*>(p);
            return getInteractionWithSuperQuad(pQuad, timeStamp, interactionHandler);
        }
    }
    return nullptr;
}

/*!
 * Computes all relevant information for the contact between this and a superquadric particle. Returns nullptr if there
 * is no (positive) overlap.
 * \param[in] p                     Superquadric particle we want to have the contact-information for
 * \param[in] timeStamp             Time stamp to be assigned to the interaction object (i.e., the current time)
 * \param[in] interactionHandler    BaseInteraction container from where the interaction is retrieved, and to which it
 *                                  is assigned (if it is a new interaction).
 * \return                          A vector of Interaction* with size 1 (if there is an interaction) or 0
 *                                  (if there is no interaction).
 */
BaseInteraction* SuperQuadricParticle::getInteractionWithSuperQuad(SuperQuadricParticle* const p, const unsigned timeStamp,
                                          InteractionHandler* const interactionHandler)
{
    BaseInteraction* const C = interactionHandler->getInteraction(p, this, timeStamp);
    const SmallVector<4> approximateContactPoint = getContactPoint(p, C);
    
    //set the new contact point:
    LabFixedCoordinates contactPoint;
    contactPoint.X = approximateContactPoint[0];
    contactPoint.Y = approximateContactPoint[1];
    contactPoint.Z = approximateContactPoint[2];
    
    //if the minimum is positive, there is no contact: return empty vector
    if (computeShape(contactPoint) > -1e-10)
    {
        return nullptr;
    }
    
    C->setContactPoint(contactPoint);
    C->setLagrangeMultiplier(approximateContactPoint[3]);
    
    SmallVector<3> gradThis = computeShapeGradientLabFixed(contactPoint).getNormalised();
    LabFixedCoordinates normal;
    normal.X = gradThis[0];
    normal.Y = gradThis[1];
    normal.Z = gradThis[2];
    C->setNormal(normal);
    
    const Mdouble alphaI = overlapFromContactPoint(contactPoint, normal);
    const Mdouble alphaJ = p->overlapFromContactPoint(contactPoint, -normal);
    C->setOverlap(alphaI + alphaJ);
    ///\todo: find correct value for 'distance'
    C->setDistance((getPosition() - p->getPosition()).getLength() - C->getOverlap());
    
    return C;
}

/*!
 * \details Computes the contact point between this and the given superquadric particle. It first gets a first guess for
 *          the contact point, and then uses that to compute the real contact-point with Newton-iterations.
 * \param[in] p the particle with which we want to know the contact point with.
 * \param[in] C The contact between this and the give superquadric. Does not contain information of this time-step yet.
 *              This contact is used to see if there was an interaction between p and this superquadric in the time-step
 *              before, and if so, that contact point can serve as an initial point for contact-detection between the
 *              particles.
 * \return      A SmallVector of length 4, with the contact point and the lagrange-multiplier needed for Newton's method
 *              for finding the contact point.
 */
SmallVector<4> SuperQuadricParticle::getContactPoint(const SuperQuadricParticle* const p, BaseInteraction* C) const
{
    //Assuming contact between two spheres
    SmallVector<4> approximateContactPoint = getInitialGuessForContact(p, C);
    if (computeContactPoint(approximateContactPoint, this, p))
    {
        return approximateContactPoint;
    }
    else
    {
        return getContactPointPlanB(p, 4);
    }
}

/*!
 * Compute the distance between the contact-point and surface of this superquadric particle with Newton-iterations.
 * This is the procedure as described by  Comp. Part. Mech. (2017) 4 : 101-118, eq 37.
 * \param[in] contactPoint The contact point between this particle with another BaseInteractable
 * \param[in] normal        The normal direction of the contact.
 * \return Distance between contact-point and surface of this particle
 */
Mdouble SuperQuadricParticle::overlapFromContactPoint(const LabFixedCoordinates& contactPoint,
        const LabFixedCoordinates& normal) const
{
    Mdouble alphaI = 0;
    LabFixedCoordinates xEdge = contactPoint + alphaI * normal;
    Mdouble dampingCoefficient = 1;
    while (std::abs(computeShape(xEdge)) > 1e-10)
    {
        SmallVector<3> gradientShape = computeShapeGradientLabFixed(xEdge);
        LabFixedCoordinates gradient;
        gradient.X = gradientShape(0);
        gradient.Y = gradientShape(1);
        gradient.Z = gradientShape(2);
        Mdouble newValue = alphaI - dampingCoefficient * computeShape(xEdge) / Vec3D::dot(gradient, normal);
        const LabFixedCoordinates xEdgeNew = contactPoint + newValue * normal;
        
        if (std::abs(computeShape(xEdgeNew)) < std::abs(computeShape(xEdge)))
        {
            alphaI = newValue;
            xEdge = xEdgeNew;
            dampingCoefficient = 1;
        }
        else
        {
            dampingCoefficient /= 2;
            if (dampingCoefficient < 1e-10)
            {
                logger(ERROR, "overlap detection does not converge");
            }
        }
    }
    return alphaI;
}

/*!
 * \details For the contact-detection, we need an initial guess where the contact will be (Newton's method only
 * converges if you start sufficiently close). If there was already a contact during the last time-step, the values of
 * the contact-point of last time-step is taken, otherwise it is taken as the middle between both particle-positions.
 * \param pQuad other superquadric for which we want to know if there is a contact
 * \param C contact between this particle and the given particle
 * \return A guess for the initial contact: contact point X,Y,Z, lagrange multiplier. Note that these are given in the
 *          LAB-FIXED coordinates.
 */
SmallVector<4> SuperQuadricParticle::getInitialGuessForContact(const SuperQuadricParticle* pQuad, BaseInteraction* const C) const
{
    SmallVector<4> approximateContactPoint;
    if (C == nullptr || C->getOverlap() < 1e-15)
    {
        //this is a new interaction
        const LabFixedCoordinates branchVector = pQuad->getPosition() - getPosition();
        //Get the square of the distance between particle i and particle j
        const Mdouble distanceSquared = Vec3D::getLengthSquared(branchVector);
        const Mdouble distance = sqrt(distanceSquared);
        const LabFixedCoordinates normal = (branchVector / distance);
        const Mdouble overlap = getSumOfInteractionRadii(pQuad) - distance;
        const LabFixedCoordinates contactPoint = (pQuad->getPosition() -
                                                  (pQuad->getInteractionRadius(this) - 0.5 * overlap) * normal);
        
        approximateContactPoint[0] = contactPoint.X;
        approximateContactPoint[1] = contactPoint.Y;
        approximateContactPoint[2] = contactPoint.Z;
        approximateContactPoint[3] = 1;
    }
    else
    {
        //this is an existing interaction
        const LabFixedCoordinates contactPoint = C->getContactPoint();
        const Mdouble multiplier = C->getLagrangeMultiplier();
        approximateContactPoint[0] = contactPoint.X;
        approximateContactPoint[1] = contactPoint.Y;
        approximateContactPoint[2] = contactPoint.Z;
        approximateContactPoint[3] = multiplier;
    }
    return approximateContactPoint;
}

/*!
 * \details This function computes the value of the shape function in the lab-fixed coordinate system. The expression
 * is provided in Section 2.3 of the article in Comp. Part. Mech. (2017) 4 : 101-118.
 * \return The value of the shape-function at the given (lab-fixed) coordinates.
 */
Mdouble SuperQuadricParticle::computeShape(const LabFixedCoordinates& labFixedCoordinates) const
{
    BodyFixedCoordinates bodyFixedCoordinates = labFixedCoordinates - this->getPosition();
    getOrientation().rotateBack(bodyFixedCoordinates);
    
    const Mdouble n1 = 2. / eps1_;
    const Mdouble n2 = 2. / eps2_;
    
    return std::pow(std::pow(std::abs(bodyFixedCoordinates.X / axes_.X), n2)
                    + std::pow(std::abs(bodyFixedCoordinates.Y / axes_.Y), n2), n1 / n2)
           + std::pow(std::abs(bodyFixedCoordinates.Z / axes_.Z), n1) - 1.0;
}

/*!
 * \details This function computes the gradient ("first derivative") of the shape function in the lab-fixed coordinate
 * system. The expressions are provided in Eq. 14 of the article in Comp. Part. Mech. (2017) 4 : 101-118.
 * \return The gradient of the shape function at the given (lab-fixed) coordinates. Note, that this is the gradient to
 * the lab-fixed coordinates, \nabla_X F(X)
 * \todo Come up with good expression for when x = y = 0 and n1 < n2
 */
SmallVector<3> SuperQuadricParticle::computeShapeGradientLabFixed(const LabFixedCoordinates& labFixedCoordinates) const
{
    
    BodyFixedCoordinates bodyFixedCoordinates = labFixedCoordinates - this->getPosition();
    getOrientation().rotateBack(bodyFixedCoordinates);
    
    const Mdouble n1 = 2. / eps1_;
    const Mdouble n2 = 2. / eps2_;
    
    const Mdouble absXoa = std::abs(bodyFixedCoordinates.X / axes_.X);
    const Mdouble absYob = std::abs(bodyFixedCoordinates.Y / axes_.Y);
    const Mdouble absZoc = std::abs(bodyFixedCoordinates.Z / axes_.Z);
    
    const Mdouble nu = std::pow(absXoa, n2) + std::pow(absYob, n2);
    
    const Mdouble help1 = std::pow(nu, n1 / n2 - 1.0);
    
    Vec3D gradientBodyFixed = {
            (n1 / axes_.X) * std::pow(absXoa, n2 - 1.0) * help1 * mathsFunc::sign(bodyFixedCoordinates.X),
            (n1 / axes_.Y) * std::pow(absYob, n2 - 1.0) * help1 * mathsFunc::sign(bodyFixedCoordinates.Y),
            (n1 / axes_.Z) * std::pow(absZoc, n1 - 1.0) * mathsFunc::sign(bodyFixedCoordinates.Z)};
    //we get the gradient in terms of the lab-fixed coordinates by \nabla_X F(X) = A \nabla_x f(x).
    getOrientation().rotate(gradientBodyFixed);
    return {gradientBodyFixed.X, gradientBodyFixed.Y, gradientBodyFixed.Z};
}

/*!
 * \details This function computes the Hessian ("second derivative") of the shape function in the body-fixed coordinate
 * system. The expressions are provided in Eq. 15 of the article in Comp. Part. Mech. (2017) 4 : 101-118.
 * \return The hessian of the shape function at the given (lab-fixed) coordinates. Note, that this is the hessian to
 * the lab-fixed coordinates, H_X (F)(X)
 * \todo Come up with good expression for when x = y = 0 and n1 < n2
 */
SmallMatrix<3, 3> SuperQuadricParticle::computeHessianLabFixed(const LabFixedCoordinates& labFixedCoordinates) const
{
    SmallMatrix<3, 3> hessian;
    BodyFixedCoordinates bodyFixedCoordinates = labFixedCoordinates - this->getPosition();
    getOrientation().rotateBack(bodyFixedCoordinates);
    
    const Mdouble n1 = 2.0 / eps1_;
    const Mdouble n2 = 2.0 / eps2_;
    
    const Mdouble absXoa = std::abs(bodyFixedCoordinates.X / axes_.X);
    const Mdouble absYob = std::abs(bodyFixedCoordinates.Y / axes_.Y);
    const Mdouble absZoc = std::abs(bodyFixedCoordinates.Z / axes_.Z);
    
    const Mdouble nu = std::pow(absXoa, n2) + std::pow(absYob, n2);
    const Mdouble help1 = std::pow(nu, n1 / n2 - 1.0);
    const Mdouble help2 = std::pow(nu, n1 / n2 - 2.0);
    
    hessian(0, 0) = n1 * (n2 - 1) / axes_.X / axes_.X * std::pow(absXoa, n2 - 2.0) * help1
                    + ((n1 - n2) * n1 / axes_.X / axes_.X) * std::pow(absXoa, 2. * n2 - 2.0) * help2;
    hessian(1, 1) = n1 * (n2 - 1) / axes_.Y / axes_.Y * std::pow(absYob, n2 - 2.0) * help1
                    + ((n1 - n2) * n1 / axes_.Y / axes_.Y) * std::pow(absYob, 2. * n2 - 2.0) * help2;
    hessian(2, 2) = n1 * (n1 - 1) / axes_.Z / axes_.Z * std::pow(absZoc, n1 - 2.0);
    hessian(1, 0) = n1 * (n1 - n2) / axes_.X / axes_.Y * std::pow(absXoa, n2 - 1.0) * std::pow(absYob, n2 - 1.0) *
                    help2 * mathsFunc::sign(bodyFixedCoordinates.X * bodyFixedCoordinates.Y);
    hessian(0, 1) = hessian(1, 0);
    SmallMatrix<3, 3> A;
    getOrientation().getRotationMatrix(A);
    return A * hessian * (A.transpose());
}

/*!
 * For the contact detection, we formulate the optimisation problem "minimise F1 + F2, s.t. F1 = F2", where F1 and F2
 * are the shape functions of particles 1 and 2. Define a Lagrange multiplier mu^2, then this is equivalent to solving
 * "Gradient1 + mu^2 Gradient2 = 0, F1 - F2 = 0". The left-hand side of that function is computed in this function. See
 * also Comp. Part. Mech. (2017) 4 : 101-118, equation 22.
 * \param[in] position  Current guess for the contact-point, where the expression above should be evaluated
 * \param[in] p1        First particle for which we are looking for a contact with.
 * \param[in] p2        Second particle for which we are looking for a contact with.
 * \return              Residual of the objective function.
 */
SmallVector<4> SuperQuadricParticle::computeResidualContactDetection(const SmallVector<4>& position,
                                                             const SuperQuadricParticle* const p1,
                                                             const SuperQuadricParticle* const p2) const
{
    LabFixedCoordinates labFixedCoordinates;
    labFixedCoordinates.X = position[0];
    labFixedCoordinates.Y = position[1];
    labFixedCoordinates.Z = position[2];
    const Mdouble lagrangeMultiplier = position[3];
    
    // First compute the gradient of the superquadric shape function and then rotate it to the
    // global frame of reference
    SmallVector<3> gradThis = p1->computeShapeGradientLabFixed(labFixedCoordinates);
    
    SmallVector<3> gradOther = p2->computeShapeGradientLabFixed(labFixedCoordinates);
    
    // Eqn. (23) in Podhlozhnyuk et al, Comp. Part. Mech. 2017.
    // where contactPoint[3] = mu
    SmallVector<4> toBecomeZero;
    for (unsigned int i = 0; i < 3; ++i)
    {
        toBecomeZero[i] = gradThis[i] + lagrangeMultiplier * lagrangeMultiplier * gradOther[i];
    }
    toBecomeZero[3] = p1->computeShape(labFixedCoordinates) - p2->computeShape(labFixedCoordinates);
    
    return toBecomeZero;
}

/*!
 * \details Compute and return the derivative of computeResidualContactDetection, both to the position
 * and the Lagrange multiplier, and evaluated at the contact point. This is done in order to use Newton's method on
 * computeResidualContactDetection, which comes from Comp. Part. Mech. (2017) 4 : 101-118.
 *
 * The result is a 4x4 matrix, with 4 parts:
 * [0,2]x[0,2] (upper-left corner): Hessian1 + mu^2 Hessian2, where Hessian is the "second derivative" of the shape
 *                                  function and mu^2 the Lagrange multiplier. This is the derivative of (22a) to x in
 *                                  Comp. Part. Mech. (2017) 4 : 101-118.
 * [0,2]x[3] (upper-right corner):  2 * mu * Gradient2, where Gradient is the gradient of the shape function and mu is
 *                                  the Lagrange multiplier. This is the derivative of (22a) to mu in Comp. Part. Mech.
 *                                  (2017) 4 : 101-118.
 * [3]x[0,2] (lower-left corner):   Gradient1 - Gradient2,  where Gradient is the gradient of the shape function.
 *                                  This is the derivative of (22b) to x in Comp. Part. Mech. (2017) 4 : 101-118.
 * [3]x[3] (lower-right corner):    0. This is the derivative of (22b) to mu in Comp. Part. Mech. (2017) 4 : 101-118.
 *
 * In order to compute these parts, first compute the first and second derivative of the shape function (the gradient
 * and hessian) in the lab-fixed coordinate system for both particles.
 * To do so, we apply the rotation matrix A on the body-fixed first and second derivative of the shape function,
 * which are computed in the functions computeShapeGradientBodyFixed() and computeHessian(). The expressions for the
 * purpose of mapping are provided in Eq. 18 of the article in Comp. Part. Mech. (2017) 4 : 101-118.
 * Then fill in the respective contributions.
 *
 * \param[in]contactPoint   The contact point where this Jacobian should be evaluated
 * \param[in] p1            First particle in the contact for which this Jacobian should be evaluated
 * \param[in] p2            Second particle in the contact for which this Jacobian should be evaluated
 * \return                  A 4x4 matrix with the Jacobian of computeResidualContactDetection.
 */
SmallMatrix<4, 4>
SuperQuadricParticle::getJacobianOfContactDetectionObjective(const SmallVector<4>& contactPoint,
                                                     const SuperQuadricParticle* const p1,
                                                     const SuperQuadricParticle* const p2) const
{
    LabFixedCoordinates labFixedCoordinates;
    labFixedCoordinates.X = contactPoint[0];
    labFixedCoordinates.Y = contactPoint[1];
    labFixedCoordinates.Z = contactPoint[2];
    const Mdouble lagrangeMultiplier = contactPoint[3];
    
    SmallVector<3> gradThis = p1->computeShapeGradientLabFixed(labFixedCoordinates);
    SmallVector<3> gradOther = p2->computeShapeGradientLabFixed(labFixedCoordinates);
    
    SmallMatrix<3, 3> hessianThis = p1->computeHessianLabFixed(labFixedCoordinates);
    SmallMatrix<3, 3> hessianOther = p2->computeHessianLabFixed(labFixedCoordinates);
    
    SmallMatrix<3, 3> upperLeftCorner = hessianThis + lagrangeMultiplier * lagrangeMultiplier * hessianOther;
    
    SmallMatrix<4, 4> jacobian;
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            jacobian(i, j) = upperLeftCorner(i, j);
        }
        jacobian(i, 3) = 2 * lagrangeMultiplier * gradOther[i];
    }
    
    for (unsigned int j = 0; j < 3; ++j)
    {
        jacobian(3, j) = gradThis[j] - gradOther[j];
    }
    
    return jacobian;
}

/*!
 * Compute the mean curvature, see Comp. Part. Mech. (2017) 4 : 101-118, eq (39)
 * \param[in] labFixedCoordinates position in lab fixed coordinate system
 * \return mean curvature of particle at labFixedCoordinates
 */
Mdouble SuperQuadricParticle::getCurvature(const LabFixedCoordinates& labFixedCoordinates) const
{
    SmallVector<3> gradientVec = computeShapeGradientLabFixed(labFixedCoordinates);
    SmallMatrix<3, 1> gradient = gradientVec;
    SmallMatrix<3, 3> hessian = computeHessianLabFixed(labFixedCoordinates);
    Mdouble helper = ((gradient.transpose() * hessian) * gradient)(0, 0);
    Mdouble gradientLength = gradientVec.length();
    Mdouble helper2 = gradientLength * gradientLength * (hessian(0, 0) + hessian(1, 1) + hessian(2, 2));
    Mdouble denominator = 2 * gradientLength * gradientLength * gradientLength;
    
    ///\todo Minus sign the wrong way around, check why
    return (helper2 - helper) / denominator;
}

/*!
 * \details Get whether or not this superquadric is in contact with the given particle: first transform the particle to
 * a superquadric if necessary, then compute the contact-point using the function getContactPoint. If the shape-function
 * at the contact point is negative, then there is a contact with the given particle, and otherwise not
 * \param[in] p The particle for which we want to know if there is a contact
 * \return      True if there is a contact, false otherwise.
 */
bool SuperQuadricParticle::isInContactWith(const BaseParticle* const p) const
{
    SmallVector<4> approximateContactPoint;

    if (p->isSphericalParticle())
    {
        SuperQuadricParticle* pQuad = new SuperQuadricParticle;
        pQuad->setAxes(p->getRadius(), p->getRadius(), p->getRadius());
        approximateContactPoint = getContactPoint(pQuad, nullptr);
        delete pQuad;
    } else {
        const SuperQuadricParticle* pQuad = static_cast<const SuperQuadricParticle*>(p);
        approximateContactPoint = getContactPoint(pQuad, nullptr);
    }

    //set the new contact point:
    LabFixedCoordinates contactPoint;
    contactPoint.X = approximateContactPoint[0];
    contactPoint.Y = approximateContactPoint[1];
    contactPoint.Z = approximateContactPoint[2];
    
    //if the minimum is positive, there is no contact: return false
    return (computeShape(contactPoint) < 0);
    
}

/*!
 * If the "normal" method of finding a contact point diverges, then we start plan B. For this, we approximate the
 * current and given particle by spheres, and move more and more towards the actual shape for the particles in small
 * increments. This way, we have a relatively good starting point for each new optimisation problem, and therefore the
 * chance that our contact-detection algorithm diverges is much smaller. This is based on equations (24) - (27) of Comp.
 * Part. Mech. (2017) 4 : 101-118.
 * \param[in] pOther    Particle we want to know if the contact-point with
 * \return              The point where the shape-functions of both particles are minimised.
 */
SmallVector<4> SuperQuadricParticle::getContactPointPlanB(const SuperQuadricParticle* const pOther, unsigned numberOfSteps) const
{
    logger(VERBOSE, "Number of iterations: %", numberOfSteps);
    // Step 1: Compute contact point for two volume equivalent spheres
    const Mdouble interactionRadiusThis = getInteractionRadius(pOther);
    const Mdouble interactionRadiusOther = pOther->getInteractionRadius(this);
    const Vec3D positionThis = getPosition();
    const Vec3D positionOther = pOther->getPosition();
    
    //analytical solution for contact point between bounding spheres
    const Vec3D initialPoint = (interactionRadiusOther * positionThis + interactionRadiusThis * positionOther) /
                               (interactionRadiusThis + interactionRadiusOther);
    
    SmallVector<4> contactPointPlanB;
    // Analytical solution
    contactPointPlanB[0] = initialPoint.X;
    contactPointPlanB[1] = initialPoint.Y;
    contactPointPlanB[2] = initialPoint.Z;
    contactPointPlanB[3] = 1.0; //Need to check.
    
    writeDebugMessageStep1(pOther, contactPointPlanB);
    
    //Step 2: Compute the deltas
    const Vec3D dAxesThis = (getAxes() - interactionRadiusThis * Vec3D(1, 1, 1)) / numberOfSteps;
    const Mdouble dn11 = (2.0 / getExponentEps1() - 2.0) / numberOfSteps;
    const Mdouble dn12 = (2.0 / getExponentEps2() - 2.0) / numberOfSteps;
    
    const Vec3D dAxesOther = (pOther->getAxes() - interactionRadiusOther * Vec3D(1, 1, 1)) / numberOfSteps;
    const Mdouble dn21 = (2.0 / pOther->getExponentEps1() - 2.0) / numberOfSteps;
    const Mdouble dn22 = (2.0 / pOther->getExponentEps2() - 2.0) / numberOfSteps;
    
    writeDebugMessageStep2(pOther, dAxesThis, dn11, dn12, dAxesOther, dn21, dn22);
    
    // Create two superquadrics with the above parameters for incrementally computing the contact point
    SuperQuadricParticle p1;
    SuperQuadricParticle p2;
    
    p1.setPosition(getPosition());
    p2.setPosition(pOther->getPosition());
    
    p1.setOrientation(getOrientation());
    p2.setOrientation(pOther->getOrientation());
    bool success = true;
    const unsigned maxIterations = 20;
    for (unsigned int counter = 1; counter <= numberOfSteps; ++counter)
    {
        
        // Step 3: Increment the shape and blockiness parameters
        const Vec3D axesThis = interactionRadiusThis * Vec3D(1, 1, 1) + counter * dAxesThis;
        const Mdouble n11 = 2.0 + counter * dn11;
        const Mdouble n12 = 2.0 + counter * dn12;
        
        Vec3D axesOther = interactionRadiusOther * Vec3D(1, 1, 1) + counter * dAxesOther;
        const Mdouble n21 = 2.0 + counter * dn21;
        const Mdouble n22 = 2.0 + counter * dn22;
        
        p1.setAxesAndExponents(axesThis, 2.0 / n11, 2.0 / n12);
        p2.setAxesAndExponents(axesOther, 2.0 / n21, 2.0 / n22);
        
        writeDebugMessageStep3(axesThis, n11, n12, axesOther, n21, n22);
        
        writeDebugMessageMiddleOfLoop(p1, p2, contactPointPlanB, counter);
        
        //compute the contact point of the new particles
        success = computeContactPoint(contactPointPlanB, &p1, &p2);
        if (!success)
        {
            if (numberOfSteps > maxIterations)
            {
                write(std::cout);
                pOther->write(std::cout);
                logger(ERROR, "Plan B fails even with more than 20 intermediate steps");
            }
            return getContactPointPlanB(pOther, 2 * numberOfSteps);
        }
    }
    logger.assert(p1.getAxes().X == getAxes().X, "In getContactPointPlanB, final particle for contact-detection not "
            "the same as original particle");
    
    logger(VERBOSE, "Plan B contact point: %", contactPointPlanB);
    
    return contactPointPlanB;
}

/*!
 * \details Function that actually performs the Newton iterations for contact-detection. We use a Newton-method with
 * adaptive damping, i.e. we start with an undamped Newton iteration, and if we "overshoot" we use a damping-factor that
 * gets halved until there is no overshoot anymore. If the damping-factor gets too small, either plan B is initiated, or
 * if this is already plan B, we give an error. After each iteration, the damping factor is set to 1 again.
 * \param[in|out] contactPoint  The contact point we are looking for. The input is an approximate contact point that is
 *                              obtained by another method, and we update this point until we reach the defined
 *                              tolerance.
 * \param[in] p1                The first particle of the contact we are looking for
 * \param[in] p2                The second particle of the contact we are looking for
 * \return                      Boolean for whether or not the contact-detection was successful.
 */
bool SuperQuadricParticle::computeContactPoint(SmallVector<4>& contactPoint, const SuperQuadricParticle* const p1,
                                       const SuperQuadricParticle* const p2) const
{
    // Damped Newton's method: (dampingFactor 1 is undamped)
    Mdouble dampingFactor = 1;
    int iteration = 1;
    const Mdouble tolerance = 1e-5;
    const unsigned maxIterations = 100;
    
    logger(VERBOSE, "func to be zero value: % \n",
           computeResidualContactDetection(contactPoint, p1, p2));
    
    while (computeResidualContactDetection(contactPoint, p1, p2).length() > tolerance && iteration < maxIterations)
    {
        logger(VERBOSE, "Iteration: %\n", iteration);
        
        SmallMatrix<4, 4> jacobian = getJacobianOfContactDetectionObjective(contactPoint, p1, p2);
        SmallVector<4> func = -computeResidualContactDetection(contactPoint, p1, p2);
        jacobian.solve(func); //note: solve is in place
        
        SmallVector<4> newPoint = contactPoint + dampingFactor * func;
        const auto residueOld = computeResidualContactDetection(contactPoint, p1, p2);
        const auto residueNew = computeResidualContactDetection(newPoint, p1, p2);
        logger(VERBOSE, "Old value: % (%) new value: % (%)",
               computeResidualContactDetection(contactPoint, p1, p2),
               computeResidualContactDetection(contactPoint, p1, p2).length(),
               computeResidualContactDetection(newPoint, p1, p2),
               computeResidualContactDetection(newPoint, p1, p2).length());
        
        logger(VERBOSE, "current contact point: %, new contact point: %\n", contactPoint, newPoint);
        logger(VERBOSE, "damping factor: %", dampingFactor);
        
        if (residueNew.length()
            < residueOld.length())
        {
            contactPoint = newPoint;
            dampingFactor = 1;
        }
        else
        {
            dampingFactor /= 2;
            
            if (dampingFactor < 1e-10)
            {
                return false;
            }
        }
        
        iteration++;
    }
    return (iteration != maxIterations);
}

void SuperQuadricParticle::writeDebugMessageMiddleOfLoop(const SuperQuadricParticle& p1, const SuperQuadricParticle& p2,
                                                 SmallVector<4>& contactPointPlanB, const unsigned int& counter) const
{
    logger(VERBOSE, "Position of particle 1 (p1): % \nPosition of particle 2 (p2): % \n",
           p1.getPosition(), p2.getPosition());
    logger(VERBOSE, "Orientation of particle 1 (p1): % \nOrientation of particle 2 (p2): %\n",
           p1.getOrientation(), p2.getOrientation());
    
    // Step 4: Calculate the contact point for the updated shape parameters
    
    logger(VERBOSE, "----------------------------------------");
    logger(VERBOSE, "STEP 4: Compute contact point");
    logger(VERBOSE, "----------------------------------------");
    logger(VERBOSE, " ");
    logger(VERBOSE, "Counter: %", counter);
    logger(VERBOSE, " ");
    
    logger(VERBOSE, "Analytical solution Contact Point: % ", contactPointPlanB);
    logger(VERBOSE, " ");
}

void
SuperQuadricParticle::writeDebugMessageStep3(const Vec3D& axesThis, const Mdouble& n11, const Mdouble& n12,
                                     const Vec3D& axesOther, const Mdouble& n21, const Mdouble& n22) const
{
    logger(VERBOSE, "-----------------------------------");
    logger(VERBOSE, "STEP 3: First increment");
    logger(VERBOSE, "-----------------------------------");
    logger(VERBOSE, " ");
    
    logger(VERBOSE, "Particle 1 x-scale after increment: %", axesThis.X);
    logger(VERBOSE, "Particle 1 y-scale after increment: %", axesThis.Y);
    logger(VERBOSE, "Particle 1 z-scale after increment: %", axesThis.Z);
    logger(VERBOSE, "Particle 1 n1 after increment: %", n11);
    logger(VERBOSE, "Particle 1 n2 after increment: %", n12);
    logger(VERBOSE, " ");
    
    logger(VERBOSE, "Particle 2 x-scale after increment: %", axesOther.X);
    logger(VERBOSE, "Particle 2 y-scale after increment: %", axesOther.Y);
    logger(VERBOSE, "Particle 2 z-scale after increment: %", axesOther.Z);
    logger(VERBOSE, "Particle 2 n1 after increment: %", n21);
    logger(VERBOSE, "Particle 2 n2 after increment: %", n22);
    logger(VERBOSE, " ");
}

void
SuperQuadricParticle::writeDebugMessageStep2(const SuperQuadricParticle* pQuad, const Vec3D& dAxesThis, const Mdouble& dn11,
                                     const Mdouble& dn12, const Vec3D& dAxesOther, const Mdouble& dn21,
                                     const Mdouble& dn22) const
{
    logger(VERBOSE, "---------------------");
    logger(VERBOSE, "STEP 2");
    logger(VERBOSE, "---------------------");
    logger(VERBOSE, " ");
    
    logger(VERBOSE, "Particle 1 x-scale: %", pQuad->getAxes().X);
    logger(VERBOSE, "Particle 1 y-scale: %", pQuad->getAxes().Y);
    logger(VERBOSE, "Particle 1 z-scale: %", pQuad->getAxes().Z);
    logger(VERBOSE, "Particle 1 n1: %", 2. / pQuad->getExponentEps1());
    logger(VERBOSE, "Particle 1 n2: %", 2. / pQuad->getExponentEps2());
    logger(VERBOSE, " ");
    
    logger(VERBOSE, "Particle 1 x-scale increment: %", dAxesThis.X);
    logger(VERBOSE, "Particle 1 y-scale increment: %", dAxesThis.Y);
    logger(VERBOSE, "Particle 1 z-scale increment: %", dAxesThis.Z);
    logger(VERBOSE, "Particle 1 n1 increment: %", dn11);
    logger(VERBOSE, "Particle 1 n2 increment: %", dn12);
    logger(VERBOSE, " ");
    
    logger(VERBOSE, "Particle 2 x-scale: %", getAxes().X);
    logger(VERBOSE, "Particle 2 y-scale: %", getAxes().Y);
    logger(VERBOSE, "Particle 2 z-scale: %", getAxes().Z);
    logger(VERBOSE, "Particle 2 n1: %", 2. / getExponentEps1());
    logger(VERBOSE, "Particle 2 n2: %", 2. / getExponentEps2());
    logger(VERBOSE, " ");
    
    logger(VERBOSE, "Particle 2 x-scale increment: %", dAxesOther.X);
    logger(VERBOSE, "Particle 2 y-scale increment: %", dAxesOther.Y);
    logger(VERBOSE, "Particle 2 z-scale increment: %", dAxesOther.Z);
    logger(VERBOSE, "Particle 2 n1 increment: %", dn21);
    logger(VERBOSE, "Particle 2 n2 increment: %", dn22);
    logger(VERBOSE, " ");
}

void SuperQuadricParticle::writeDebugMessageStep1(const SuperQuadricParticle* pQuad, const SmallVector<4>& contactPointPlanB) const
{
    logger(VERBOSE, "---------------------");
    logger(VERBOSE, "STEP 1");
    logger(VERBOSE, "---------------------\n");
    
    logger(VERBOSE, "Position of particle 1: % ", getPosition());
    logger(VERBOSE, "Position of particle 2 (pQuad): % ", pQuad->getPosition());
    logger(VERBOSE, " ");
    
    logger(VERBOSE, "Radius of particle 1: % ", getInteractionRadius(pQuad));
    logger(VERBOSE, "Radius of particle 2 (pQuad): % \n", pQuad->getInteractionRadius(this));
    
    logger(VERBOSE, "Orientation of particle 1: % ", getOrientation());
    logger(VERBOSE, "Orientation of particle 2 (pQuad): % ", pQuad->getOrientation());
    logger(VERBOSE, " ");
    
    logger(VERBOSE, "Particle 1 axes: %", getAxes());
    logger(VERBOSE, "Particle 2 axes (pQuad): % \n", pQuad->getAxes());
    
    logger(VERBOSE, "Analytical solution for two equivalent spheres in contact: % \n", contactPointPlanB);
}

Mdouble SuperQuadricParticle::getInteractionRadius(const BaseParticle* particle) const
{
    const auto mixedSpecies = getSpecies()->getHandler()->getMixedObject(getSpecies(),particle->getSpecies());
    return getRadius() + 0.5 * mixedSpecies->getInteractionDistance();
}

/// Computes the particle's (inverse) mass and inertia.
void SuperQuadricParticle::computeMass(const ParticleSpecies& s) {
    if (isFixed()) return;
    if (getParticleDimensions()==3) {
        Mdouble volume = getVolume();

        Mdouble help1 = mathsFunc::beta(1.5 * eps2_, 0.5 * eps2_);
        Mdouble help2 = mathsFunc::beta(0.5 * eps1_, 2.0 * eps1_ + 1.0);
        Mdouble help3 = mathsFunc::beta(0.5 * eps2_, 0.5 * eps2_ + 1);
        Mdouble help4 = mathsFunc::beta(1.5 * eps1_, eps1_ + 1.0);

        invMass_ = 1.0 / (volume * s.getDensity());
        invInertia_.XX = 1.0 / (s.getDensity() * (0.5 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) *
                                   (axes_.Y * axes_.Y * help1 * help2
                                    + 4.0 * axes_.Z * axes_.Z * help3 * help4));
        invInertia_.XY = 0.0;
        invInertia_.XZ = 0.0;
        invInertia_.YY = 1.0 / (s.getDensity() * (0.5 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) *
                                   (axes_.X * axes_.X * help1 * help2
                                    + 4.0 * axes_.Z * axes_.Z * help3 * help4));
        invInertia_.YZ = 0.0;
        invInertia_.ZZ = 1.0 / (s.getDensity() * (0.5 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) *
                                   (axes_.X * axes_.X + axes_.Y * axes_.Y) * help1 * help2);
    } else {
        logger(ERROR, "[SuperQuadricParticle::computeMass] SuperQuadricParticle cannot be two-dimensional (yet)");
    }
};
