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

#include <string>
#include <Mercury3D.h>
#include "ClumpParticle.h"
#include "ParticleHandler.h"
#include "DPMBase.h"

#define N_ITER 3 // Parameter of time integration
// Lower values give better performance, higher-better precision.
// (Passing self-test requires at least 2, PFC4.0 uses 4)


ClumpParticle::ClumpParticle()
{
    setRadius(1.0);
    nPebble_ = 0;

    clumpMass_ = 1.0;
    viscousDamping_ = 0.0;
    pebblePos_ = std::vector<Vec3D>(0);

    pebbleRadius_ = std::vector<Mdouble>(0);
    pebbleParticles_ = std::vector<ClumpParticle*>(0);

    principalDirections_ = Matrix3D(1, 0, 0, 0, 1, 0, 0, 0, 1);
    initPrincipalDirections_ = Matrix3D(1, 0, 0, 0, 1, 0, 0, 0, 1);

    clumpInertia_ = MatrixSymmetric3D(1, 0, 0, 1, 0, 1);
    clumpInitInertia_ = MatrixSymmetric3D(1, 0, 0, 1, 0, 1);
    rotationMatrix_ = Matrix3D(1, 0, 0, 0, 1, 0, 0, 0, 1);
    invInertia_= clumpInertia_.inverse();

    isPebble_ = false; //Assign false by default
    isClump_ = false; //Assign false by default
    clumpParticle = nullptr;


    DzhanibekovParticle_ = false;
    VerticallyOriented_ = false;

    logger(DEBUG, "Clump() created");
}

ClumpParticle::ClumpParticle(const ClumpParticle& p): NonSphericalParticle(p)
{
    nPebble_ = p.nPebble_;

    clumpMass_ = p.clumpMass_;
    viscousDamping_ = p.viscousDamping_;
    pebblePos_ = p.pebblePos_;
    pebbleRadius_ = p.pebbleRadius_;
    pebbleParticles_ = p.pebbleParticles_;
    principalDirections_ = p.principalDirections_;
    initPrincipalDirections_ = p.initPrincipalDirections_;
    clumpInertia_ = p.clumpInertia_;
    clumpInitInertia_ = p.clumpInitInertia_;
    invInertia_= clumpInertia_.inverse();
    DzhanibekovParticle_ = p.DzhanibekovParticle_;
    VerticallyOriented_ = p.VerticallyOriented_;

    for (int iPebble = 1; iPebble <= nPebble_; iPebble++) pebbleParticles_[iPebble - 1] = nullptr;

    // Pebble attributes
    isPebble_ = p.isPebble_;
    clumpParticle = p.clumpParticle;
    // Clump attributes
    isClump_ = p.isClump_;
}

ClumpParticle::~ClumpParticle()
= default;


ClumpParticle* ClumpParticle::copy() const
{
    return new ClumpParticle(*this);
}


std::string ClumpParticle::getName() const
{
    return "ClumpParticle";
}

void ClumpParticle::read(std::istream& is)
{
    BaseParticle::read(is);
    std::string dummy;
    is >> dummy >> nPebble_;
}

int ClumpParticle::NPebble() const
{
    return nPebble_;
}

void ClumpParticle::setClump()
{
    isClump_ = true;
    isPebble_ = false;


}

//Function to store pebble information
void ClumpParticle::addPebble(Vec3D position, Mdouble radius)
{
    nPebble_++; //Counter of pebbles
    pebblePos_.push_back(position); //Store pebble positions
    pebbleRadius_.push_back(radius); //Store pebble radius
    pebbleParticles_.push_back(nullptr);//Store null pointer per pebble.
//    isClump_ = true;
}

// V
void ClumpParticle::setPrincipalDirections(Matrix3D directions)
{
    principalDirections_ = directions;
}

// V
void ClumpParticle::setInitPrincipalDirections(Matrix3D directions)
{
    initPrincipalDirections_ = directions;
}

// Rotate basis vectors around given rotation vector by |rotation|
void ClumpParticle::rotatePrincipalDirections(Vec3D rotation)
{
    Mdouble tol = 10e-9;
    Mdouble angle = rotation.getLength();
    if (angle < tol) return;
    Mdouble c1, c2, c3, theta, dist;
    Vec3D e3 = rotation; e3.normalise();
    Vec3D e2 = Vec3D::cross(e3, getPrincipalDirections_e1());
    if (e2.getLength() < 0.1) e2 = Vec3D::cross(e3, getPrincipalDirections_e2());
    if (e2.getLength() < 0.1) e2 = Vec3D::cross(e3, getPrincipalDirections_e3());
    e2.normalise();
    Vec3D e1 = Vec3D::cross(e2, e3);
    Vec3D pa;
    // Rotate principal axis 1
    pa = getPrincipalDirections_e1();
    c1 = Vec3D::dot(pa, e1);
    c2 = Vec3D::dot(pa, e2);
    c3 = Vec3D::dot(pa, e3);

    dist = sqrt(c1 * c1 + c2 * c2);
    theta = atan2(c2, c1) + angle;
    setPrincipalDirections_e1(dist * cos(theta) * e1 + dist*sin(theta) * e2 + c3 * e3);
    // Rotate principal axis 2
    pa = getPrincipalDirections_e2();
    c1 = Vec3D::dot(pa, e1);
    c2 = Vec3D::dot(pa, e2);
    c3 = Vec3D::dot(pa, e3);
    dist = sqrt(c1 * c1 + c2 * c2);
    theta = atan2(c2, c1) + angle;
    setPrincipalDirections_e2(dist * cos(theta) * e1 + dist*sin(theta) * e2 + c3 * e3);
    // Rotate principal axis 3
    pa = getPrincipalDirections_e3();
    c1 = Vec3D::dot(pa, e1);
    c2 = Vec3D::dot(pa, e2);
    c3 = Vec3D::dot(pa, e3);
    dist = sqrt(c1 * c1 + c2 * c2);
    theta = atan2(c2, c1) + angle;
    setPrincipalDirections_e3(dist * cos(theta) * e1 + dist*sin(theta) * e2 + c3 * e3);
}

//After adding the objects to the handler, update clump quantities.
void ClumpParticle::actionsAfterAddObject()
{
    //Only attribute features to clump particle.
    if (isPebble()) return;
    if (isClump())
    {
        ClumpParticle p0;   // Instance for the pebbles
        p0.setSpecies(getSpecies());
        p0.setClump(this);

        //Go through the number of pebbles
        for (int iPebble = 1; iPebble <= NPebble(); iPebble++)
        {
            //Set the radius of pebbles:
            p0.setRadius(getPebbleRadius()[iPebble - 1]);
            //Store the address of the pebble:
            setPebble(iPebble - 1, getHandler()->copyAndAddObject(p0));
        }

        //Update the pebble positions, velocities.
        updatePebblesVelPos();
    }
}

void ClumpParticle::setInitInertia(MatrixSymmetric3D inertia)
{
    clumpInitInertia_ = inertia;
    clumpInertia_ = inertia;
    invInertia_= clumpInertia_.inverse();
}

// Rotate TOI
void ClumpParticle::rotateTensorOfInertia()
{
   // Initial and current principal directions
   Vec3D e10 = getInitPrincipalDirections_e1();
   Vec3D e20 = getInitPrincipalDirections_e2();
   Vec3D e30 = getInitPrincipalDirections_e3();

   Vec3D e1 = getPrincipalDirections_e1();
   Vec3D e2 = getPrincipalDirections_e2();
   Vec3D e3 = getPrincipalDirections_e3();

   // Rotation matrix from initial and current principal directions
   Matrix3D Q(Vec3D::dot(e10, e1), Vec3D::dot(e10, e2), Vec3D::dot(e10, e3),
              Vec3D::dot(e20, e1), Vec3D::dot(e20, e2), Vec3D::dot(e20, e3),
              Vec3D::dot(e30, e1), Vec3D::dot(e30, e2), Vec3D::dot(e30, e3));

   Matrix3D Qt = transpose(Q);

  MatrixSymmetric3D inertia = MtoS(Q * (StoM(clumpInitInertia_) * Qt));
  //  inertia = clumpInitInertia_;  // uncomment to turn off rotation of toi
  clumpInertia_ = inertia;
  invInertia_= clumpInertia_.inverse();
  rotationMatrix_ = Q;

}


// Rotate pebbles
void ClumpParticle::updatePebblesVelPos()
{
    BaseParticle* pPebble;
    Vec3D position = getPosition();
    Vec3D angularVelocity = getAngularVelocity();
    Quaternion orientation = getOrientation();

    Vec3D e1 = getPrincipalDirections_e1();
    Vec3D e2 = getPrincipalDirections_e2();
    Vec3D e3 = getPrincipalDirections_e3();

    Vec3D velocityDueToRotation;

    for (int iPebble = 1; iPebble <= nPebble_; iPebble++)
    {
        pPebble = pebbleParticles_[iPebble - 1];
        // pPebble->setMass(clumpMass_);
        pPebble->invMass_ = 1. / clumpMass_;
        pPebble->setAngularVelocity(angularVelocity);
        pPebble->setOrientation(orientation);
        pPebble->setPosition(position + e1 * pebblePos_[iPebble - 1].X + e2 * pebblePos_[iPebble - 1].Y + e3 * pebblePos_[iPebble - 1].Z);


        velocityDueToRotation = Vec3D::cross(angularVelocity, pPebble->getPosition() - position);

        pPebble->setVelocity(getVelocity());
        pPebble->addVelocity(velocityDueToRotation);
    }
}

/*!
 * \details First step of Velocity Verlet integration (see also
 * http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet).
 * \param[in] time          current time
 * \param[in] timeStep      current time step
 */
void ClumpParticle::integrateBeforeForceComputation(double time, double timeStep)
{
    if (isPebble()) return;
    if (getInvMass() == 0.0)
    {
        BaseInteractable::integrateBeforeForceComputation(time, timeStep);
    }
    else
    {
#ifdef MERCURY_USE_MPI
        //For periodic particles in parallel the previous position is required
        setPreviousPosition(getPosition());
#endif
        accelerate((getForce() - viscousDamping_ * getVelocity()) * getInvMass() * 0.5 * timeStep); // V(t+0.5dt)
        const Vec3D displacement = getVelocity() * timeStep;
        move(displacement); // X(t+dt)

        DPMBase* const dpm = getHandler()->getDPMBase();
        if (!dpm->getHGridUpdateEachTimeStep())
        {
            dpm->hGridUpdateMove(this, displacement.getLengthSquared());
        }

        // PFC4 style acceleration of Clumps
        angularAccelerateClumpIterative(timeStep); //W(t+0.5dt)

        //apply to rotation quaternion q: q = normalise(q + \tilde{C}\omega*timeStep) (see Wouter's notes)
        rotate(getAngularVelocity() * timeStep);

        // Rotate
        rotatePrincipalDirections(getAngularVelocity() * timeStep);

        // Update pebble nodes
        updatePebblesVelPos();

        // Update TOI
        rotateTensorOfInertia();

    }
}

/*!
 * \details Second step of Velocity Verlet integration (see also
 * http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet).
 * \param[in] time      current time
 * \param[in] timeStep  current time step
 */
void ClumpParticle::integrateAfterForceComputation(double time, double timeStep)
{
    if (isPebble()) return;
    if (getInvMass() == 0.0)
    {
        BaseInteractable::integrateAfterForceComputation(time, timeStep);
    }
    else
    {
        // Translational acceleration
        accelerate((getForce() - viscousDamping_ * getVelocity()) * getInvMass() * 0.5 * timeStep);

        // PFC4 style angular acceleration of Clumps
        angularAccelerateClumpIterative(timeStep);

        // Rotate pebbles
        updatePebblesVelPos();

        // Updates extra quantities
        updateExtraQuantities();
    }
}

// Clump-specific time integration algorithm (PFC4.0 style)
void ClumpParticle::angularAccelerateClumpIterative(double timeStep)
{
    Mdouble Ixx_ = getInertia().XX;
    Mdouble Iyy_ = getInertia().YY;
    Mdouble Izz_ = getInertia().ZZ;
    Mdouble Ixy_ = getInertia().XY;
    Mdouble Ixz_ = getInertia().XZ;
    Mdouble Iyz_ = getInertia().YZ;

    Vec3D M = getTorque() - viscousDamping_ * getAngularVelocity();

    Vec3D angularVelocity_0 = getAngularVelocity();
    Vec3D angularVelocity_n = angularVelocity_0;
    for (int n = 0; n<N_ITER; n++) //
    {
        Mdouble wx = angularVelocity_n.X;
        Mdouble wy = angularVelocity_n.Y;
        Mdouble wz = angularVelocity_n.Z;

        Vec3D W = Vec3D(wy*wz*(Izz_-Iyy_) - wz*wz*Iyz_ + wy*wy*Iyz_ + wx*wy*Ixz_ - wz*wx*Ixy_,
                        wz*wx*(Ixx_-Izz_) - wx*wx*Ixz_ + wz*wz*Ixz_ + wy*wz*Ixy_ - wx*wy*Iyz_,
                        wx*wy*(Iyy_-Ixx_) - wy*wy*Ixy_ + wx*wx*Ixy_ + wz*wx*Iyz_ - wy*wz*Ixz_);

        Vec3D angularAcceleration_n = invInertia_ * (M-W);
        angularAcceleration_ = angularAcceleration_n;
        angularAcceleration_ = angularAcceleration_n;
        angularVelocity_n = angularVelocity_0 + 0.5 * timeStep * angularAcceleration_n;
    }

    setAngularVelocity(angularVelocity_n);
}

// V
void ClumpParticle::computeMass(const ParticleSpecies &s)
{
    if (isFixed()) return;
    if (isPebble()) return;

    if (isClump())
    {
        invMass_ = 1.0 / clumpMass_;
        invInertia_= clumpInertia_.inverse();

    }
}


void ClumpParticle::updateExtraQuantities()
{
    Mdouble ANG_TOL = 0.1; // 5.7 degrees - tolerance to misalignment
    Mdouble ACC_TOL = 0.1; // Tolerance for angular acceleration magnitude

    Mdouble TOL = 10e-8;   // External force tolerance
    ClumpParticle* pPebble;
    Vec3D n3 = getPrincipalDirections_e3();
    Vec3D v = Vec3D(0,0,1);

    // Check for vertical alignment
    if (acos(Vec3D::dot(n3, v))<ANG_TOL)
    {
        setVerticallyOriented(true);
        for (int iPebble = 1; iPebble <= nPebble_; iPebble++)
        {
            pPebble = pebbleParticles_[iPebble - 1];
            pPebble->setVerticallyOriented(true);
        }
    }
    else {
        setVerticallyOriented(false);
        for (int iPebble = 1; iPebble <= nPebble_; iPebble++)
        {
            pPebble = pebbleParticles_[iPebble - 1];
            pPebble->setVerticallyOriented(false);
        }
    }

    // Check for Dzhanibekov States

    Vec3D w  = getAngularVelocity()/getAngularVelocity().getLength();
    Vec3D n2 = getPrincipalDirections_e2();
    Mdouble acc = angularAcceleration_.getLength();

    if ((acos(Vec3D::dot(n2, w))<ANG_TOL) || (acc<TOL))
    {
        setDzhanibekovParticle(true);
        for (int iPebble = 1; iPebble <= nPebble_; iPebble++)
        {
            pPebble = pebbleParticles_[iPebble - 1];
            pPebble->setDzhanibekovParticle(true);
        }
    }

    // Any contact force/torques break D state
    if ( (getForce().getLength() > TOL)||(getTorque().getLength() > TOL) )
    {
        setDzhanibekovParticle(false);
        for (int iPebble = 1; iPebble <= nPebble_; iPebble++)
        {
            pPebble = pebbleParticles_[iPebble - 1];
            pPebble->setDzhanibekovParticle(false);
        }
    }

    return;
}