//Copyright (c) 2013-2022, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef MultiParticle_H
#define MultiParticle_H

#include "Mercury3D.h"
#include "Math/SmallVector.h"
#include "BaseParticle.h"
#include "BaseInteractable.h"
#include "NonSphericalParticle.h"
#include "Species/ParticleSpecies.h"
#include "DPMBase.h"
/*!
 * \class MultiParticle
 * \brief
 */
class MultiParticle final: public NonSphericalParticle
{
public:
    /*!
     * \brief Basic Particle constructor, creates a particle at (0,0,0) with radius, mass and inertia equal to 1
     */
    MultiParticle();

    /*!
     * \brief Copy constructor, which accepts as input a reference to a Superquadric.
     * It creates a copy of this Particle and all it's information.
     * Usually it is better to use the copy() function for polymorphism.
     */
    MultiParticle(const MultiParticle& p);
    
    /*!
    * \brief Destructor, needs to be implemented and checked to see if it is the largest or smallest particle currently
    * in its particleHandler
    */

    ~MultiParticle() override;

    MultiParticle* copy() const override;

    void write(std::ostream& os) const override;

    void read(std::istream& is) override;

    std::string getName() const override;

    bool isSphericalParticle() const override
    {
        return true;
    }

    void computeMass(const ParticleSpecies& s) override;

    void setInitInertia(MatrixSymmetric3D inertia);

    void rotateTensorOfInertia();

    void addSlave(Vec3D position, Mdouble radius);

    void setMaster();

    void setPrincipalDirections(Matrix3D directions);

    Vec3D getPrincipalDirections_e1() const
    {
        return Vec3D(principalDirections.XX, principalDirections.YX, principalDirections.ZX);
    }
    Vec3D getPrincipalDirections_e2() const
    {
        return Vec3D(principalDirections.XY, principalDirections.YY, principalDirections.ZY);
    }
    Vec3D getPrincipalDirections_e3() const
    {
        return Vec3D(principalDirections.XZ, principalDirections.YZ, principalDirections.ZZ);
    }

    // Methods to obtain initial principal directions
    Vec3D getInitPrincipalDirections_e1() const
    {
        return Vec3D(initPrincipalDirections.XX, initPrincipalDirections.YX, initPrincipalDirections.ZX);
    }
    Vec3D getInitPrincipalDirections_e2() const
    {
        return Vec3D(initPrincipalDirections.XY, initPrincipalDirections.YY, initPrincipalDirections.ZY);
    }
    Vec3D getInitPrincipalDirections_e3() const
    {
        return Vec3D(initPrincipalDirections.XZ, initPrincipalDirections.YZ, initPrincipalDirections.ZZ);
    }


    int NSlave() const;

    void actionsAfterAddObject() override;

    void updateSlavesVelPos();

    void updateSlavesVel();

    void integrateBeforeForceComputation(double time, double timeStep) override;

    void integrateAfterForceComputation(double time, double timeStep) override;

    void angularAccelerateMasterIterative(double timeStep);

    void rotatePrincipalDirections(Vec3D rotation);

    // Principle direction
    void setPrincipalDirections_e1(Vec3D e)
    {
        principalDirections.XX = e.X;
        principalDirections.YX = e.Y;
        principalDirections.ZX = e.Z;
    }
    void setPrincipalDirections_e2(Vec3D e)
    {
        principalDirections.XY = e.X;
        principalDirections.YY = e.Y;
        principalDirections.ZY = e.Z;
    }
    void setPrincipalDirections_e3(Vec3D e)
    {
        principalDirections.XZ = e.X;
        principalDirections.YZ = e.Y;
        principalDirections.ZZ = e.Z;
    }

    std::vector<Mdouble> getSlaveRadius() const {
        return slaveRadius;
    }

    void setSlave(int kSlave, MultiParticle* pSlave) {
        slaveParticles[kSlave] = pSlave;
    }

    //ToDo: This function is used in ParticleHandler
    void setMaster(MultiParticle* master) {
        isMaster = false;
        isSlave = true;
        masterParticle = master;
    }

    void setMassMultiparticle(Mdouble mass)
    {
        massMultiparticle = mass;
    }

    void setDamping(Mdouble damp)
    {
        viscousDamping = damp;
    }

private:

    int nSlave;

    Mdouble massMultiparticle;
    Mdouble viscousDamping;
    MatrixSymmetric3D inertiaMultiparticle;
    MatrixSymmetric3D initInertiaMultiparticle;

    std::vector<Vec3D> slavePos;
    std::vector<Mdouble> slaveRadius;

    Matrix3D principalDirections;
    Matrix3D initPrincipalDirections;
    std::vector<MultiParticle*> slaveParticles;

    //Helper functions
    MatrixSymmetric3D MtoS( Matrix3D M){ return MatrixSymmetric3D(M.XX, M.XY, M.XZ, M.YY, M.YZ, M.ZZ);}
    Matrix3D StoM( MatrixSymmetric3D M){ return Matrix3D(M.XX, M.XY, M.XZ, M.XY, M.YY, M.YZ, M.XZ, M.YZ, M.ZZ);}
    Matrix3D transpose(Matrix3D M){ return Matrix3D(M.XX, M.YX, M.ZX, M.XY, M.YY, M.ZY, M.XZ, M.YZ, M.ZZ);}
};

#endif
