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


    void write(std::ostream& os) const override
    {
        BaseParticle::write(os);
        os << " DzhanibekovParticle " << DzhanibekovParticle_;
        os << " VerticallyOriented " << VerticallyOriented_;

    }

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

    // Principal direction
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
        invMass_ = 1 / massMultiparticle;
    }

    void setDamping(Mdouble damp)
    {
        viscousDamping = damp;
    }

    Mdouble getKineticEnergy() const override{
        Mdouble res = 0;
        if (isMaster) {
            Vec3D v = getVelocity();
            res = 0.5 * massMultiparticle * ( v.X * v.X +  v.Y * v.Y + v.Z * v.Z );
        }
        return res;
    }

    Mdouble getRotationalEnergy() const override{
        Mdouble res = 0;
        if (isMaster) {
            Vec3D nn = getAngularVelocity();
            Mdouble nl = nn.getLength();
            Mdouble tol = 1e-10;
            if (nl > tol) {
                nn /= nl;
                Mdouble ww = getAngularVelocity().getLengthSquared();
                Mdouble II = Vec3D::dot(nn, (getInertia() * nn));
                res = 0.5 * II * ww;
            }
        }
        return res;

    }

    std::vector <Vec3D> getSlavePositions(){
        std::vector <Vec3D> globalPos;
        Vec3D e1 = getPrincipalDirections_e1();
        Vec3D e2 = getPrincipalDirections_e2();
        Vec3D e3 = getPrincipalDirections_e3();
        for (int i = 1; i<=NSlave(); i++){
        globalPos.push_back(getPosition() + e1*slavePos[i-1].X + e2*slavePos[i-1].Y + e3*slavePos[i-1].Z);
        }
        return globalPos;
    }

    std::vector <Mdouble> getSlaveRadii(){ return slaveRadius; }

    // Methods setting and getting some extra Boolean properties

    bool getDzhanibekovParticle()
    {
    	return DzhanibekovParticle_;
    }
    
    bool getVerticallyOriented()
    {
    	return VerticallyOriented_;
    }
    
    void setDzhanibekovParticle( bool d)
    {
    	DzhanibekovParticle_ = d;
    }
    
    void setVerticallyOriented( bool d)
    {
    	VerticallyOriented_ = d;
    }

    unsigned getNumberOfFieldsVTK() const override
    {
        return 2;
    }

    std::string getTypeVTK(unsigned i) const override
    {
        return "Int8";
    }


    std::string getNameVTK(unsigned i) const override
    {
        if (i==0)
            return "DzhanibekovParticle";
        else
            return "VerticallyOriented";
    }


    std::vector<Mdouble> getFieldVTK(unsigned i) const override
    {
        if (i==0)
            return std::vector<Mdouble>(1, DzhanibekovParticle_);
        else
            return std::vector<Mdouble>(1, VerticallyOriented_);
    }

    void updateExtraQuantities();
		

private:

    int nSlave;
    
    bool DzhanibekovParticle_; // This property is needed to quantify Dzhanibekov gas properties
    bool VerticallyOriented_;  // This property is useful for mechnical stability simulations (Gomboc, Dominos)
    Vec3D angularAcceleration_;

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
