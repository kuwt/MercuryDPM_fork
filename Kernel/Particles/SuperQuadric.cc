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

#include "SuperQuadric.h"
#include "InteractionHandler.h"
#include "ParticleHandler.h"

/*!
 * \details calls the default constructor of BaseParticle, creates an SuperEllipsoid at (0,0,0) assuming it is
 * sphere.
 */
SuperQuadric::SuperQuadric()
{
    axes_ = Vec3D(1.0, 1.0, 1.0);
    eps1_ = 1.0;
    eps2_ = 1.0;
    
    //setInertia();
    
    logger(DEBUG, "SuperQuadric::SuperQuadric() finished");
}

/*!
 * \details Constructor that copies most of the properties of the given particle.
 *          Please note that not everything is copied, for example the position 
 *          in the HGrid is not determined yet by the end of this constructor. 
 *          It also does not copy the interactions and the pointer to the handler
 *          that handles this particle. Use with care.
 * \param[in,out] p  Reference to the SuperQuad this one should become a copy of.
 */
SuperQuadric::SuperQuadric(const SuperQuadric& p)
        : BaseParticle(p)
{
    axes_ = p.axes_;
    eps1_ = p.eps1_;
    eps2_ = p.eps2_;
}

/*!
 * \details Destructor. It asks the ParticleHandler to check if this was the 
 *          smallest or largest particle and adjust itself accordingly.
 */
SuperQuadric::~SuperQuadric()
{
    if (getHandler() != nullptr)
    {
        getHandler()->checkExtremaOnDelete(this);
    }
    logger(DEBUG, "SuperQuadric::SuperQuadricric() of particle % finished.", getId());
    
}

/*!
 * \details Copy method. Uses copy constructor to create a copy on the heap. 
 *          Useful for polymorphism.
 * \return pointer to the particle's copy
 */
SuperQuadric* SuperQuadric::copy() const
{
    return new SuperQuadric(*this);
}

/*!
 * \details SuperQuad print method, which accepts an os std::ostream as
 *          input. It prints human readable SuperQuad information to the
 *          std::ostream.
 * \param[in,out] os    stream to which the info is written
 */
void SuperQuadric::write(std::ostream& os) const
{
    BaseParticle::write(os);
    os << " axes " << axes_
       << " exp1 " << eps1_
       << " exp2 " << eps2_;
    
}

/*!
 * \details Returns the name of the object; in this case 'SuperQuad'.
 * \return The object name.
 */
std::string SuperQuadric::getName() const
{
    return "SuperQuadric";
}


/*!
 * \details Particle read function, which reads the axes_ and both epsilons from the given input-stream.
 * \param[in,out] is    input stream with particle properties.
 */
void SuperQuadric::read(std::istream& is)
{
    BaseParticle::read(is);
    std::string dummy;
    is >> dummy >> axes_ >> dummy >> eps1_ >> dummy >> eps2_;
    logger.assert_always(eps1_ < 2 + 1e-10, "epsilon1 should be at most 2");
    logger.assert_always(eps2_ < 2 + 1e-10, "epsilon2 should be at most 2");
}

// Set and Get Functions

void SuperQuadric::setAxesAndExponents(const Mdouble& a1, const Mdouble& a2, const Mdouble& a3, const Mdouble& eps1,
                                    const Mdouble& eps2)
{
    axes_.X = a1;
    axes_.Y = a2;
    axes_.Z = a3;
    eps1_ = eps1;
    eps2_ = eps2;
    logger.assert_always(eps1_ < 2 + 1e-10, "epsilon1 should be at most 2");
    logger.assert_always(eps2_ < 2 + 1e-10, "epsilon2 should be at most 2");
}

void SuperQuadric::setAxesAndExponents(const Vec3D& axes, const Mdouble& eps1, const Mdouble& eps2)
{
    axes_ = axes;
    eps1_ = eps1;
    eps2_ = eps2;
    logger.assert_always(eps1_ < 2 + 1e-10, "epsilon1 should be at most 2");
    logger.assert_always(eps2_ < 2 + 1e-10, "epsilon2 should be at most 2");
}

void SuperQuadric::setAxes(const Vec3D& axes)
{
    axes_ = axes;
}

void SuperQuadric::setAxes(const Mdouble& a1, const Mdouble& a2, const Mdouble& a3)
{
    axes_.X = a1;
    axes_.Y = a2;
    axes_.Z = a3;
}

void SuperQuadric::setExponents(const Mdouble& eps1, const Mdouble& eps2)
{
    eps1_ = eps1;
    eps2_ = eps2;
    logger.assert_always(eps1_ < 2 + 1e-10, "epsilon1 should be at most 2");
    logger.assert_always(eps2_ < 2 + 1e-10, "epsilon2 should be at most 2");
}

Vec3D SuperQuadric::getAxes()
{
    return axes_;
}

Mdouble SuperQuadric::getExponentEps1()
{
    logger.assert_always(eps1_ < 2 + 1e-10, "epsilon1 should be at most 2");
    return eps1_;
}

Mdouble SuperQuadric::getExponentEps2()
{
    logger.assert_always(eps2_ < 2 + 1e-10, "epsilon2 should be at most 2");
    return eps2_;
}

/*!
 * \details Returns the volume of the SuperEllipsoid, which is calculated using its principal axis and
 * the exponents. The analytical expressions are taken from Chapter 2 of the book "Segmentation and Recovery of Superquadrics"
 * by Jaklic et al. However, the beta functions that are part of these expressions are approximations, see ExtendedMath.cc
 * \return The actual volume of this SuperEllipsoid.
 */
Mdouble SuperQuadric::getVolume() const
{
    if (getHandler() == nullptr)
    {
        logger(ERROR, "[SuperQuadric::getVolume] no particle handler specified");
        return 0;
    }
    
    return (2.0 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) * mathsFunc::beta(0.5 * eps1_ + 1.0, eps1_) *
           mathsFunc::beta(0.5 * eps2_, 0.5 * eps2_);
    
}

/*!
 * \details This function computes the principal moments of inertia for the superellipsoids. Again, the analytical expressions
 * are taken from Chapter 2 of the book "Segmentation and Recovery of Superquadrics" by Jaklic et al.
 */
void SuperQuadric::setInertia()
{
    MatrixSymmetric3D inertia;
    
    inertia.XX = getSpecies()->getDensity() * (0.5 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) *
                 (axes_.Y * axes_.Y * mathsFunc::beta(1.5 * eps2_, 0.5 * eps2_) * mathsFunc::beta(0.5 * eps1_, 2.0 * eps1_ + 1.0)
                  + 4.0 * axes_.Z * axes_.Z * mathsFunc::beta(0.5 * eps2_, 0.5 * eps2_ + 1) *
                    mathsFunc::beta(1.5 * eps1_, eps1_ + 1.0));//0.2 * mass * (axes_.Y * axes_.Y + axes_.Z * axes_.Z);
    inertia.YY = getSpecies()->getDensity() * (0.5 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) *
                 (axes_.X * axes_.X * mathsFunc::beta(1.5 * eps2_, 0.5 * eps2_) * mathsFunc::beta(0.5 * eps1_, 2.0 * eps1_ + 1.0)
                  + 4.0 * axes_.Z * axes_.Z * mathsFunc::beta(0.5 * eps2_, 0.5 * eps2_ + 1) *
                    mathsFunc::beta(1.5 * eps1_, eps1_ + 1.0));//0.2 * mass * (axes_.X * axes_.X + axes_.Z * axes_.Z);
    inertia.ZZ = getSpecies()->getDensity() * (0.5 * axes_.X * axes_.Y * axes_.Z * eps1_ * eps2_) *
                 (axes_.X * axes_.X + axes_.Y * axes_.Y) * mathsFunc::beta(1.5 * eps2_, 0.5 * eps2_) *
                 mathsFunc::beta(0.5 * eps1_, 2.0 * eps1_ + 1.0); //0.2 * mass * (axes_.Y * axes_.Y + axes_.X * axes_.X);
    
    BaseParticle::setInertia(inertia);
}

/*!
 * \details This function determines the minimum bounding radius of the bounding sphere. Bounding sphere is a bounding
 * volume technique where a simple sphere encapsulates a more complex shaped object. Additionally, we also for other
 * surrounding effects such as liquid films or electromagnetic fields. For spheres and ellipsoids the mininmum bounding
 * radius is max(a1,a2,a3), where (a1,a2,a3) are the principal axes. For generalised superellipsoidic shapes, we determine
 * the minimum radius by utilising the solution presented in Eq. 12 of the article by Podlozhynuk et al. in Comp. Part. Mech.
 * (2017). They basically solve an optimisation problem using Lagrange multipliers. However note that in their solutions n1 = 2.0/eps1_
 * and n2 = 2.0/eps2_
 * \return The radius of the bounding sphere + 0.5 *interactionDistance (in case of e.g. liquid films around the
 *          particle)
 */
Mdouble SuperQuadric::getInteractionRadius() const
{
    // For spheres and ellipsoids
    if (mathsFunc::isEqual(eps1_, 1, 1.e-16) && mathsFunc::isEqual(eps2_, 1, 1e-16))
    {
        const Mdouble boundingRadius = std::max(std::max(axes_.X, axes_.Y), axes_.Z);
        return boundingRadius + getSpecies()->getInteractionDistance() * 0.5;
    }
        // For cylinders
    else if (mathsFunc::isEqual(eps2_, 1.0, 1e-16))
    {
        ///\todo automatically let a > b or implement the case b >= a
        ///\todo check validity for a = b
        logger.assert(axes_.X >= axes_.Y, "Need a >= b for bounding sphere computation. Please re-arrange your axes.");
        const Mdouble beta = std::pow(axes_.Z * axes_.Z / (axes_.X * axes_.X), 1.0 / (2.0 / eps1_ - 2.0));
        const Mdouble xTilde = std::pow(1.0 + std::pow(beta, 2.0 / eps1_), -eps1_ / 2);
        const Mdouble boundingRadius = std::sqrt(mathsFunc::square(axes_.X * xTilde)
                                                 + mathsFunc::square(beta * axes_.Z * xTilde));
        return boundingRadius + getSpecies()->getInteractionDistance() * 0.5;
    }
    else if (mathsFunc::isEqual(eps1_, 1.0, 1e-16))
    {
        const Mdouble alpha = std::pow(axes_.Y / axes_.X, 2.0 / (2.0 / eps2_ - 2.0));
        const Mdouble help1 = std::pow(alpha, 2.0 / eps2_);
        const Mdouble gamma = std::pow(1.0 + help1, eps2_ / eps1_ - 1.0);
        logger.assert(gamma * axes_.Z * axes_.Z < axes_.X * axes_.X, "Please choose your axes such that a >= c");
        const Mdouble xTilde = (std::pow((1 + help1), -eps2_ / 2.0));
        const Mdouble boundingRadius = std::sqrt(mathsFunc::square(axes_.X * xTilde)
                                                 + mathsFunc::square(alpha * axes_.Y * xTilde));
        return boundingRadius + getSpecies()->getInteractionDistance() * 0.5;
    }
    else  // For other shapes
    {
        const Mdouble alpha = std::pow(axes_.Y / axes_.X, 2.0 / (2.0 / eps2_ - 2.0));
        const Mdouble help1 = std::pow(alpha, 2.0 / eps2_);
        const Mdouble gamma = std::pow(1.0 + help1, eps2_ / eps1_ - 1.0);
        const Mdouble beta = std::pow(gamma * axes_.Z * axes_.Z / (axes_.X * axes_.X), 1.0 / (2.0 / eps1_ - 2.0));
        const Mdouble xTilde = std::pow(std::pow(1 + help1, eps2_ / eps1_) + std::pow(beta, 2.0 / eps1_),
                                        -eps1_ / 2.0);
        const Mdouble boundingRadius = std::sqrt(mathsFunc::square(axes_.X * xTilde)
                                                 + mathsFunc::square(alpha * axes_.Y * xTilde)
                                                 + mathsFunc::square(beta * axes_.Z * xTilde));
        return boundingRadius + getSpecies()->getInteractionDistance() * 0.5;
    }
    
}

/*!
 * Overwrites BaseInteractable::getInteractionWith. First checks if the bounding radii overlap, and if so, changes the
 * given particle into a superquadric and calls getInteractionWithSuperQuad to obtain the relevant contact data.
 * @param p
 * @param timeStamp
 * @param interactionHandler
 * @return
 */
std::vector<BaseInteraction*> SuperQuadric::getInteractionWith(
        BaseParticle* const p, const unsigned timeStamp,
        InteractionHandler* const interactionHandler)
{
    //get the normal (from P away from the contact)
    const LabFixedCoordinates branchVector = p->getPosition() - getPosition();
    //Get the square of the distance between particle i and particle j
    const Mdouble distanceSquared = Vec3D::getLengthSquared(branchVector);
    const Mdouble sumOfInteractionRadii = p->getInteractionRadius() + getInteractionRadius();
    std::vector<BaseInteraction*> contacts;
    if (distanceSquared < (sumOfInteractionRadii * sumOfInteractionRadii))
    {
        //make a superquadric out of the particle.
        SuperQuadric* pQuad = dynamic_cast<SuperQuadric*>(p);
        //if the dynamic casting of the particle pointer p is a null pointer, it implies that the p is a sphere.
        bool fromSphere = (pQuad == nullptr);
        // For the sake of simplicity, we make a superquad out of it.
        if (fromSphere)
        {
            pQuad = new SuperQuadric;
            pQuad->setAxes(p->getRadius(), p->getRadius(), p->getRadius());
        }
        
        contacts = getInteractionWithSuperQuad(pQuad, timeStamp, interactionHandler);
        
        if (fromSphere)
        {
            delete pQuad;
        }
        
    }
    return contacts;
}

/*!
 * Computes all relevant information for the contact between this and a superquadric particle. Returns nullptr if there
 * is no (positive) overlap.
 * @param pQuad
 * @param timeStamp
 * @param interactionHandler
 * @return
 */
std::vector<BaseInteraction*> SuperQuadric::getInteractionWithSuperQuad(SuperQuadric* const pQuad, const unsigned timeStamp,
                                                        InteractionHandler* const interactionHandler)
{
    BaseInteraction* const C = interactionHandler->getInteraction(pQuad, this, timeStamp);
    const SmallVector<4> approximateContactPoint = getContactPoint(pQuad, C);
    
    //set the new contact point:
    LabFixedCoordinates contactPoint;
    contactPoint.X = approximateContactPoint[0];
    contactPoint.Y = approximateContactPoint[1];
    contactPoint.Z = approximateContactPoint[2];
    
    //if the minimum is positive, there is no contact: return empty vector
    if (computeShape(contactPoint) > -1e-10)
    {
        return {};
    }
    
    C->setContactPoint(contactPoint);
    C->setLagrangeMultiplier(approximateContactPoint[3]);
    
    SmallVector<3> gradThis = computeShapeGradientLocal(contactPoint).getNormalised();
    getOrientation().rotate(gradThis);
    SmallVector<3> gradOther = pQuad->computeShapeGradientLocal(contactPoint).getNormalised();
    pQuad->getOrientation().rotate(gradOther);
    LabFixedCoordinates normal;
    normal.X = gradThis[0];
    normal.Y = gradThis[1];
    normal.Z = gradThis[2];
    C->setNormal(normal);
    
    const Mdouble alphaI = computeOverlapAlpha(contactPoint, normal);
    const Mdouble alphaJ = pQuad->computeOverlapAlpha(contactPoint, -normal);
    C->setOverlap(alphaI + alphaJ);
    
    return {C};
}

/*!
 * Computes the contact point between this and the given superquadric particle
 * @param pQuad the particle with which we want to know the contact point with.
 * @param C The contact between this and the give superquadric. Does not contain information of this time-step yet.
 * @return the contact point, with lagrange-multiplier
 */
SmallVector<4> SuperQuadric::getContactPoint(const SuperQuadric* const pQuad, BaseInteraction* C) const
{
    SmallVector<4> approximateContactPoint = getInitialGuessForContact(pQuad, C);
    
    //Damped Newton's method: (dampingFactor 1 is undamped)
    Mdouble dampingFactor = 1;
    while (functionThatShouldBecomeZeroForContactDetection(approximateContactPoint, pQuad).length() > 1e-5)
    {
        
        SmallMatrix<4, 4> jacobian = computeJacobian(approximateContactPoint, pQuad);
        SmallVector<4> func = -functionThatShouldBecomeZeroForContactDetection(approximateContactPoint, pQuad);
        jacobian.solve(func); //note: solve is in place
        SmallVector<4> newPoint = approximateContactPoint + dampingFactor * func;
        logger(DEBUG, "Old value: % new value: %",
               functionThatShouldBecomeZeroForContactDetection(approximateContactPoint, pQuad),
               functionThatShouldBecomeZeroForContactDetection(newPoint, pQuad));
        
        logger(DEBUG, "current contact point: %, new contact point: %", approximateContactPoint, newPoint);
        
        if (functionThatShouldBecomeZeroForContactDetection(newPoint, pQuad).length()
            < functionThatShouldBecomeZeroForContactDetection(approximateContactPoint, pQuad).length())
        {
            approximateContactPoint = newPoint;
            dampingFactor = 1;
        }
        else
        {
            dampingFactor /= 2;
            if (dampingFactor < 1e-10)
            {
                logger(ERROR, "Contact detection diverges");
            }
        }
    }
    return approximateContactPoint;
}

/*!
 * Compute the distance between the contact-point and surface of this superquadric particle.
 * @param contactPoint
 * @param normal
 * @return Distance between contact-point and surface
 */
Mdouble SuperQuadric::computeOverlapAlpha(const LabFixedCoordinates& contactPoint, const LabFixedCoordinates& normal) const
{
    Mdouble alphaI = 0;
    LabFixedCoordinates xEdge = contactPoint + alphaI * normal;
    Mdouble dampingCoefficient = 1;
    while (std::abs(computeShape(xEdge)) > 1e-10)
    {
        SmallVector<3> gradientShape = computeShapeGradientLocal(xEdge);
        getOrientation().rotate(gradientShape);
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
 * For the contact-detection, we need an initial guess where the contact will be (Newton's method only converges if you
 * start sufficiently close). If there was already a contact during the last time-step, the values of the last time-step
 * \param pQuad other superquadric for which we want to know if there is a contact
 * \param C contact between this particle and the given particle
 * \return A guess for the initial contact: contact point X,Y,Z, lagrange multiplier
 */
SmallVector<4> SuperQuadric::getInitialGuessForContact(const SuperQuadric* pQuad, BaseInteraction* const C) const
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
        const Mdouble overlap = (pQuad->getInteractionRadius() + getInteractionRadius() - distance);
        const LabFixedCoordinates contactPoint = (pQuad->getPosition() -
                                                (pQuad->getInteractionRadius() - 0.5 * overlap) * normal);
        
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
 * \details This function computes the value of the shape function in the local coordinate system. The expression
 * is provided in Section 2.3 of the article in Comp. Part. Mech. (2017) 4 : 101-118.
 * @return
 */
Mdouble SuperQuadric::computeShape(const LabFixedCoordinates& xGlobal) const
{
    BodyFixedCoordinates xLocal = xGlobal - this->getPosition();
    getOrientation().rotateBack(xLocal);
    
    const Mdouble n1 = 2. / eps1_;
    const Mdouble n2 = 2. / eps2_;
    
    return std::pow(std::pow(std::abs(xLocal.X / axes_.X), n2) + std::pow(std::abs(xLocal.Y / axes_.Y), n2), n1 / n2)
           + std::pow(std::abs(xLocal.Z / axes_.Z), n1) - 1.0;
}

/*!
 * \details This function computes the first derivative of the shape function in the local coordinate system. The expressions
 * are provided in Eq. 14 of the article in Comp. Part. Mech. (2017) 4 : 101-118.
 * @return
 * \todo test this
 */
SmallVector<3> SuperQuadric::computeShapeGradientLocal(const LabFixedCoordinates& xGlobal) const
{
    
    BodyFixedCoordinates xLocal = xGlobal - this->getPosition();
    getOrientation().rotateBack(xLocal);
    
    const Mdouble n1 = 2. / eps1_;
    const Mdouble n2 = 2. / eps2_;
    
    const Mdouble absXoa = std::abs(xLocal.X / axes_.X);
    const Mdouble absYob = std::abs(xLocal.Y / axes_.Y);
    
    const Mdouble nu = std::pow(absXoa, n2) + std::pow(absYob, n2);
    
    const Mdouble help1 = std::pow(nu, n1 / n2 - 1.0);
    
    return {(n1 / axes_.X) * std::pow(absXoa, n2 - 1.0) * help1 * mathsFunc::sign(xLocal.X),
            (n1 / axes_.Y) * std::pow(absYob, n2 - 1.0) * help1 * mathsFunc::sign(xLocal.Y),
            (n1 / axes_.Z) * std::pow(std::abs(xLocal.Z / axes_.Z), n1 - 1.0) * mathsFunc::sign(xLocal.Z)};
}

/*!
 * \details This function computes the first derivative of the shape function in the local coordinate system. The expressions
 * are provided in Eq. 15 of the article in Comp. Part. Mech. (2017) 4 : 101-118.
 * @return
 */
SmallMatrix<3, 3> SuperQuadric::computeHessian(const LabFixedCoordinates& xGlobal) const
{
    SmallMatrix<3, 3> hessian;
    if (mathsFunc::isEqual(eps1_, 1, 1e-10) && mathsFunc::isEqual(eps2_, 1, 1e-10))
    {
        hessian(0, 0) = 2. / axes_.X / axes_.X;
        hessian(1, 1) = 2. / axes_.Y / axes_.Y;
        hessian(2, 2) = 2. / axes_.Z / axes_.Z;
    }
        ///\todo implement cases where eps1_ = 1 xor eps2_1 = 1
        ///\todo test
    else
    {
        BodyFixedCoordinates xLocal = xGlobal - this->getPosition();
        getOrientation().rotateBack(xLocal);
        
        const Mdouble n1 = 2. / eps1_;
        const Mdouble n2 = 2. / eps2_;
        
        const Mdouble absXoa = std::abs(xLocal.X / axes_.X);
        const Mdouble absYob = std::abs(xLocal.Y / axes_.Y);
        const Mdouble absZoc = std::abs(xLocal.Z / axes_.Z);
        
        const Mdouble nu = std::pow(absXoa, n2) + std::pow(absYob, n2);
        const Mdouble help1 = std::pow(nu, n1 / n2 - 1.0);
        const Mdouble help2 = std::pow(nu, n1 / n2 - 2.0);
        
        hessian(0, 0) = n1 * (n2 - 1) / axes_.X / axes_.X * std::pow(absXoa, n2 - 2.0) * help1
                        + ((n1 - n2) * n1 / axes_.X / axes_.X) * std::pow(absXoa, 2. * n2 - 2.0) * help2;
        hessian(1, 1) = n1 * (n2 - 1) / axes_.Y / axes_.Y * std::pow(absYob, n2 - 2.0) * help1
                        + ((n1 - n2) * n1 / axes_.Y / axes_.Y) * std::pow(absYob, 2. * n2 - 2.0) * help2;
        hessian(2, 2) = n1 * (n1 - 1) / axes_.Z / axes_.Z * std::pow(absZoc, n1 - 2.0);
        hessian(1, 2) = n1 * (n1 - n2) / axes_.X / axes_.Y * std::pow(absXoa, n2 - 1.0) * std::pow(absYob, n2 - 1.0) *
                        help2 * mathsFunc::sign(xLocal.X * xLocal.Y);
        if (mathsFunc::isEqual(absXoa, 0, 1e-10))
        {
            hessian(0, 0) = 0;
            hessian(1, 2) = 0;
        }
        if (mathsFunc::isEqual(absYob, 0, 1e-10))
        {
            hessian(1, 1) = 0;
            hessian(1, 2) = 0;
        }
        if (mathsFunc::isEqual(absZoc, 0, 1e-10))
        {
            hessian(2, 2) = 0;
        }
        hessian(2, 1) = hessian(1, 2);
    }
    
    return hessian;
}

SmallVector<4> SuperQuadric::functionThatShouldBecomeZeroForContactDetection(const SmallVector<4>& approximateContact,
                                                                          const SuperQuadric* const pQuad) const
{
    LabFixedCoordinates xGlobal;
    xGlobal.X = approximateContact[0];
    xGlobal.Y = approximateContact[1];
    xGlobal.Z = approximateContact[2];
    SmallVector<3> gradThis = computeShapeGradientLocal(xGlobal);
    getOrientation().rotate(gradThis);
    
    SmallVector<3> gradOther = pQuad->computeShapeGradientLocal(xGlobal);
    pQuad->getOrientation().rotate(gradOther);
    
    SmallVector<4> toBecomeZero;
    for (unsigned int i = 0; i < 3; ++i)
    {
        toBecomeZero[i] = gradThis[i] + approximateContact[3] * approximateContact[3] * gradOther[i];
    }
    toBecomeZero[3] = computeShape(xGlobal) - pQuad->computeShape(xGlobal);
    return toBecomeZero;
}

/*!
 * \details This function computes the first and second derivative of the shape function in the global coordinate
 * system. To do so, we apply the rotation matrix A on the local first and second derivative of the shape function, which are computed
 * in the functions computeShapeGradientLocal() and computeHessian(). The expressions for the purpose of mapping are provided in
 * Eq. 18 of the article in Comp. Part. Mech. (2017) 4 : 101-118.
 * @return
 */
SmallMatrix<4, 4>
SuperQuadric::computeJacobian(const SmallVector<4>& approximateContact, const SuperQuadric* const pQuad) const
{
    LabFixedCoordinates xGlobal;
    xGlobal.X = approximateContact[0];
    xGlobal.Y = approximateContact[1];
    xGlobal.Z = approximateContact[2];
    const Mdouble lagrangeMultiplier = approximateContact[3];
    
    SmallVector<3> gradThis = computeShapeGradientLocal(xGlobal);
    getOrientation().rotate(gradThis);
    SmallVector<3> gradOther = pQuad->computeShapeGradientLocal(xGlobal);
    pQuad->getOrientation().rotate(gradOther);
    
    SmallMatrix<3, 3> hessianThis = computeHessian(xGlobal);
    getOrientation().rotateTensor(hessianThis);
    SmallMatrix<3, 3> hessianOther = pQuad->computeHessian(xGlobal);
    pQuad->getOrientation().rotateTensor(hessianOther);
    
    SmallMatrix<3, 3> upperLeftCorner = hessianThis + lagrangeMultiplier * lagrangeMultiplier * hessianOther;
    
    SmallMatrix<4, 4> jacobian;
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            jacobian(i, j) = upperLeftCorner(i, j);
        }
        jacobian(i, 3) = 2 * approximateContact[3] * gradOther[i];
    }
    
    for (unsigned int j = 0; j < 3; ++j)
    {
        jacobian(3, j) = gradThis[j] - gradOther[j];
    }
    
    return jacobian;
}

/*!
 * Compute the mean curvature, see Podlozhyuk et al. (2017)
 * @param xGlobal position in global coordinate system
 * @return mean curvature of particle at xGlobal
 */
Mdouble SuperQuadric::getCurvature(const LabFixedCoordinates& xGlobal) const
{
    SmallVector<3> gradientVec = computeShapeGradientLocal(xGlobal);
    getOrientation().rotate(gradientVec);
    SmallMatrix<3, 1> gradient = gradientVec;
    SmallMatrix<3, 3> hessian = computeHessian(xGlobal);
    getOrientation().rotateTensor(hessian);
    Mdouble helper = ((gradient.transpose() * hessian) * gradient)(0, 0);
    Mdouble gradientLength = gradientVec.length();
    Mdouble helper2 = gradientLength * gradientLength * (hessian(0, 0) + hessian(1, 1) + hessian(2, 2));
    Mdouble denominator = 2 * gradientLength * gradientLength * gradientLength;
    
    ///\todo Minus sign the wrong way around, check why
    return (helper2 - helper) / denominator;
}

bool SuperQuadric::isInContactWith(const BaseParticle* const p) const
{
    //make a superquadric out of the particle.
    SuperQuadric* pQuad = dynamic_cast<SuperQuadric*>(p->copy());
    //if the dynamic casting of the particle pointer p is a null pointer, it implies that the p is a sphere.
    bool fromSphere = (pQuad == nullptr);
    // For the sake of simplicity, we make a superquad out of it.
    if (fromSphere)
    {
        pQuad = new SuperQuadric;
        pQuad->setAxes(p->getRadius(), p->getRadius(), p->getRadius());
    }
    
    const SmallVector<4> approximateContactPoint = getContactPoint(pQuad, nullptr);
    
    if (fromSphere)
    {
        delete pQuad;
    }
    
    //set the new contact point:
    LabFixedCoordinates contactPoint;
    contactPoint.X = approximateContactPoint[0];
    contactPoint.Y = approximateContactPoint[1];
    contactPoint.Z = approximateContactPoint[2];
    
    //if the minimum is positive, there is no contact: return nullptr
    return (computeShape(contactPoint) < 0);
    
}
