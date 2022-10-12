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

#include <string>
#include <Mercury3D.h>
#include "MultiParticle.h"
#include "ParticleHandler.h"
#include "DPMBase.h"

MultiParticle::MultiParticle()
{
    setRadius(1.0);
    nSlave = 0;
    // slaveMass =  std::vector<Mdouble>(0);
    // slaveVolume =  std::vector<Mdouble>(0);
    massMultiparticle = 1.0;
    viscousDamping = 0.0;
    slavePos = std::vector<Vec3D>(0);

    slaveRadius = std::vector<Mdouble>(0);
    slaveParticles = std::vector<MultiParticle*>(0);

    principalDirections = Matrix3D(1, 0, 0, 0, 1, 0, 0, 0, 1);
    initPrincipalDirections = Matrix3D(1, 0, 0, 0, 1, 0, 0, 0, 1);

    inertiaMultiparticle = MatrixSymmetric3D(1, 0, 0, 1, 0, 1);
    initInertiaMultiparticle = MatrixSymmetric3D(1, 0, 0, 1, 0, 1);
    invInertia_= inertiaMultiparticle.inverse();

    //++++++++++++++
    isSlave = false; //Assign false by default
    isMaster = false; //Assign false by default
    masterParticle = nullptr;

    //++++++++++++++
    logger(DEBUG, "Multiparticle::MultiParticle() finished");
}

MultiParticle::MultiParticle(const MultiParticle& p): NonSphericalParticle(p)
{
    nSlave = p.nSlave;
    // slaveMass = p.slaveMass;
    // slaveVolume = p.slaveVolume;
    
    massMultiparticle = p.massMultiparticle;
    viscousDamping = p.viscousDamping;
    slavePos = p.slavePos;
    slaveRadius = p.slaveRadius;
    slaveParticles = p.slaveParticles;
    principalDirections = p.principalDirections;
    initPrincipalDirections = p.initPrincipalDirections;

    inertiaMultiparticle = p.inertiaMultiparticle;
    initInertiaMultiparticle = p.initInertiaMultiparticle;
    invInertia_= inertiaMultiparticle.inverse();



    for (int iSlave = 1; iSlave <= nSlave; iSlave++) slaveParticles[iSlave - 1] = nullptr;

    //++++++++++++
    // Slave attributes
    isSlave = p.isSlave;
    masterParticle = p.masterParticle;
    // Master attributes
    isMaster = p.isMaster;
    //++++++++++++
}

MultiParticle::~MultiParticle()
= default;


MultiParticle* MultiParticle::copy() const
{
    return new MultiParticle(*this);
}

void MultiParticle::write(std::ostream& os) const
{
    BaseParticle::write(os);
    os << " nSlaves " << nSlave;
}

std::string MultiParticle::getName() const
{
    return "MultiParticle";
}

void MultiParticle::read(std::istream& is)
{
    BaseParticle::read(is);
    std::string dummy;
    is >> dummy >> nSlave;
}

int MultiParticle::NSlave() const
{
    return nSlave;
}

void MultiParticle::setMaster()
{
    isMaster = true;

}

//Function to store slave information
void MultiParticle::addSlave(Vec3D position, Mdouble radius)
{
    nSlave++; //Counter of slaves
    slavePos.push_back(position); //Store slave positions
    slaveRadius.push_back(radius); //Store slave radius
    slaveParticles.push_back(nullptr);//Store null pointer per slave.
//    isMaster = true;
}

// V
void MultiParticle::setPrincipalDirections(Matrix3D directions)
{
    principalDirections = directions;
}

// V
void MultiParticle::rotatePrincipalDirections(Vec3D rotation)
{
    Mdouble tol = 10e-10;
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

//After adding the objects to the handler, compute multiparticle quantities.
void MultiParticle::actionsAfterAddObject()
{
    //Only attribute features to master particle.
    if (IsSlave()) return;
    if (IsMaster())
    {
        MultiParticle p0;   // Instance for the pebbles
        p0.setSpecies(getSpecies());
        p0.setMaster(this);

        //Go through the number of slaves
        for (int iSlave = 1; iSlave <= NSlave(); iSlave++)
        {
            //Set the radius of slaves:
            p0.setRadius(getSlaveRadius()[iSlave - 1]);
            //Store the address of the slave:
            setSlave(iSlave - 1, getHandler()->copyAndAddObject(p0));
        }

        //Update the slave positions, velocities.
        updateSlavesVelPos();
    }
}

void MultiParticle::setInitInertia(MatrixSymmetric3D inertia)
{
    initInertiaMultiparticle = inertia;
    inertiaMultiparticle = inertia;
    invInertia_= inertiaMultiparticle.inverse();
}

void MultiParticle::rotateTensorOfInertia()
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

  MatrixSymmetric3D inertia = MtoS(Q * (StoM(initInertiaMultiparticle) * Qt));
  //  inertia = initInertiaMultiparticle;  // uncomment to turn off rotation of toi
  inertiaMultiparticle = inertia;
  invInertia_= inertiaMultiparticle.inverse();
}



void MultiParticle::updateSlavesVelPos()
{
    BaseParticle* pSlave;
    Vec3D position = getPosition();
    Vec3D angularVelocity = getAngularVelocity();
    Quaternion orientation = getOrientation();

    Vec3D e1 = getPrincipalDirections_e1();
    Vec3D e2 = getPrincipalDirections_e2();
    Vec3D e3 = getPrincipalDirections_e3();

    Vec3D velocityDueToRotation;

    for (int iSlave = 1; iSlave <= nSlave; iSlave++)
    {
        pSlave = slaveParticles[iSlave-1];
        // pSlave->setMass(massMultiparticle);
        pSlave->invMass_ = 1./massMultiparticle;
        pSlave->setAngularVelocity(angularVelocity);
        pSlave->setOrientation(orientation);
        pSlave->setPosition(position + e1*slavePos[iSlave-1].X + e2*slavePos[iSlave-1].Y + e3*slavePos[iSlave-1].Z);
        velocityDueToRotation = Vec3D::cross(angularVelocity, slavePos[iSlave-1] - position);
        pSlave->addVelocity(velocityDueToRotation);
    }
}

void MultiParticle::updateSlavesVel()
{
    BaseParticle* pSlave;
    Vec3D position = getPosition();
    Vec3D angularVelocity = getAngularVelocity();
    Quaternion orientation = getOrientation();

    Vec3D e1 = getPrincipalDirections_e1();
    Vec3D e2 = getPrincipalDirections_e2();
    Vec3D e3 = getPrincipalDirections_e3();

    Vec3D velocityDueToRotation;

    for (int iSlave = 1; iSlave <= nSlave; iSlave++)
    {
        pSlave = slaveParticles[iSlave-1];
        // pSlave->setMass(massMultiparticle);
        pSlave->invMass_ = 1./massMultiparticle;
        pSlave->setAngularVelocity(angularVelocity);
        velocityDueToRotation = Vec3D::cross(angularVelocity, slavePos[iSlave-1] - position);
        pSlave->addVelocity(velocityDueToRotation);
    }
}



/*!
 * \details First step of Velocity Verlet integration (see also
 * http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet).
 * \param[in] time          current time
 * \param[in] timeStep      current time step
 */
void MultiParticle::integrateBeforeForceComputation(double time, double timeStep)
{
    if (IsSlave()) return;
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
        accelerate((getForce() - viscousDamping*getVelocity()) * getInvMass() * 0.5 * timeStep); // V(t+0.5dt)
        const Vec3D displacement = getVelocity() * timeStep;
        move(displacement); // X(t+dt)

        DPMBase* const dpm = getHandler()->getDPMBase();
        if (!dpm->getHGridUpdateEachTimeStep())
        {
            dpm->hGridUpdateMove(this, displacement.getLengthSquared());
        }

        // PFC4 style acceleration of clumps
        angularAccelerateMasterIterative(0.5 * timeStep); //W(t+0.5dt)

        //apply to rotation quaternion q: q = normalise(q + \tilde{C}\omega*timeStep) (see Wouter's notes)
        rotate(getAngularVelocity() * timeStep);

        // Rotatestd::cout<<"w_"<<n<<" = "<<angularVelocity_n<<std::endl; the principal axis according to the angle of rotation
        rotatePrincipalDirections(getAngularVelocity() * timeStep);

        // Update slave nodes
        updateSlavesVelPos();
        rotateTensorOfInertia();

    }
}

/*!
 * \details Second step of Velocity Verlet integration (see also
 * http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet).
 * \param[in] time      current time
 * \param[in] timeStep  current time step
 */
void MultiParticle::integrateAfterForceComputation(double time, double timeStep)
{
    if (IsSlave()) return;
    if (getInvMass() == 0.0)
    {
        BaseInteractable::integrateAfterForceComputation(time, timeStep);
    }
    else
    {
        accelerate((getForce() - viscousDamping*getVelocity()) * getInvMass() * 0.5 * timeStep);
        // PFC4 style acceleration of clumps
        angularAccelerateMasterIterative(0.5 * timeStep);
        updateSlavesVelPos();
    }
}

// V
void MultiParticle::angularAccelerateMasterIterative(double timeStep)
{
    Mdouble Ixx_ = getInertia().XX;
    Mdouble Iyy_ = getInertia().YY;
    Mdouble Izz_ = getInertia().ZZ;
    Mdouble Ixy_ = getInertia().XY;
    Mdouble Ixz_ = getInertia().XZ;
    Mdouble Iyz_ = getInertia().YZ;

    Vec3D M = getTorque() - viscousDamping * getAngularVelocity();

    Vec3D angularVelocity_0 = getAngularVelocity();
    Vec3D angularVelocity_n = angularVelocity_0;
    for (int n = 0; n<1; n++) // 4 iterations needed for good precision (PFC4.0 manual)
    {
        Mdouble wx = angularVelocity_n.X;
        Mdouble wy = angularVelocity_n.Y;
        Mdouble wz = angularVelocity_n.Z;

        Vec3D W = Vec3D(wy*wz*(Izz_-Iyy_) - wz*wz*Iyz_ + wy*wy*Iyz_ + wx*wy*Ixz_ - wz*wx*Ixy_,
                        wz*wx*(Ixx_-Izz_) - wx*wx*Ixz_ + wz*wz*Ixz_ + wy*wz*Ixy_ - wx*wy*Iyz_,
                        wx*wy*(Iyy_-Ixx_) - wy*wy*Ixy_ + wx*wx*Ixy_ + wz*wx*Iyz_ - wy*wz*Ixz_);

        Vec3D angularAcceleration_n = invInertia_ * (M-W);
        angularVelocity_n = angularVelocity_0 + 0.5 * timeStep * angularAcceleration_n;
    }
    setAngularVelocity(angularVelocity_n);
}

// V
void MultiParticle::computeMass(const ParticleSpecies &s)
{
    if (isFixed()) return;
    if (IsSlave()) return;

    if (IsMaster())
    {
        invMass_ = 1.0/massMultiparticle;
        invInertia_= inertiaMultiparticle.inverse();
    }
}


