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
 * \class Clump Particle provides the functional of rigid clumps in mercuryDPM
 * \brief
 */
class ClumpParticle final: public NonSphericalParticle
{
public:
    /*!
     * \brief Basic Particle constructor, creates a particle at (0,0,0) with radius, mass and inertia equal to 1
     */
    ClumpParticle();

    /*!
    * \brief Basic Particle constructor (copy-based)
    */
    ClumpParticle(const ClumpParticle& p);
    
    /*!
    * \brief Destructor, needs to be implemented and checked to see if it is the largest or smallest particle currently
    * in its particleHandler
    */

    ~ClumpParticle() override;

    ClumpParticle* copy() const override;


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

    void addPebble(Vec3D position, Mdouble radius);

    void setClump();

    void setPrincipalDirections(Matrix3D directions);

    void setInitPrincipalDirections(Matrix3D directions);

    Vec3D getPrincipalDirections_e1() const
    {

        return Vec3D(principalDirections_.XX, principalDirections_.YX, principalDirections_.ZX);
    }
    Vec3D getPrincipalDirections_e2() const
    {
        return Vec3D(principalDirections_.XY, principalDirections_.YY, principalDirections_.ZY);
    }
    Vec3D getPrincipalDirections_e3() const
    {
        return Vec3D(principalDirections_.XZ, principalDirections_.YZ, principalDirections_.ZZ);
    }

    // Methods to obtain initial principal directions
    Vec3D getInitPrincipalDirections_e1() const
    {
        return Vec3D(initPrincipalDirections_.XX, initPrincipalDirections_.YX, initPrincipalDirections_.ZX);
    }
    Vec3D getInitPrincipalDirections_e2() const
    {
        return Vec3D(initPrincipalDirections_.XY, initPrincipalDirections_.YY, initPrincipalDirections_.ZY);
    }
    Vec3D getInitPrincipalDirections_e3() const
    {
        return Vec3D(initPrincipalDirections_.XZ, initPrincipalDirections_.YZ, initPrincipalDirections_.ZZ);
    }


    int NPebble() const; // Number of pebbles (for a clump particle)

    void actionsAfterAddObject() override; // The function that updates clump quantities after adding a pebble

    void updatePebblesVelPos();

    void integrateBeforeForceComputation(double time, double timeStep) override;

    void integrateAfterForceComputation(double time, double timeStep) override;

    void angularAccelerateClumpIterative(double timeStep); // Specific method of time integration for a clump particle

    void rotatePrincipalDirections(Vec3D rotation);

    // Principal direction
    void setPrincipalDirections_e1(Vec3D e)
    {
        principalDirections_.XX = e.X;
        principalDirections_.YX = e.Y;
        principalDirections_.ZX = e.Z;
    }
    void setPrincipalDirections_e2(Vec3D e)
    {
        principalDirections_.XY = e.X;
        principalDirections_.YY = e.Y;
        principalDirections_.ZY = e.Z;
    }
    void setPrincipalDirections_e3(Vec3D e)
    {
        principalDirections_.XZ = e.X;
        principalDirections_.YZ = e.Y;
        principalDirections_.ZZ = e.Z;
    }

    std::vector<Mdouble> getPebbleRadius() const {
        return pebbleRadius_;
    }

    // add pointer to pebble pointers list
    void setPebble(int kPebble, ClumpParticle* pPebble) {
        pebbleParticles_[kPebble] = pPebble;
    }

    // Sets the particle to be a pebble of a given clump
    void setClump(ClumpParticle* master) {
        isClump_ = false;
        isPebble_ = true;
        clumpParticle = master;
    }

    void setClumpMass(Mdouble mass)
    {
        clumpMass_ = mass;
        invMass_ = 1 / clumpMass_;
    }

    // Extra viscous damping on a clump
    void setDamping(Mdouble damp)
    {
        viscousDamping_ = damp;
    }

    Mdouble getKineticEnergy() const override{
        Mdouble res = 0;
        if (isClump_) {
            Vec3D v = getVelocity();
            res = 0.5 * clumpMass_ * (v.X * v.X + v.Y * v.Y + v.Z * v.Z );
        }
        return res;
    }

    Mdouble getRotationalEnergy() const override{
        Mdouble res = 0;
        if (isClump_) {
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

    std::vector <Vec3D> getPebblePositions(){
        std::vector <Vec3D> globalPos;
        Vec3D e1 = getPrincipalDirections_e1();
        Vec3D e2 = getPrincipalDirections_e2();
        Vec3D e3 = getPrincipalDirections_e3();
        for (int i = 1; i <= NPebble(); i++){
        globalPos.push_back(getPosition() + e1 * pebblePos_[i - 1].X + e2 * pebblePos_[i - 1].Y + e3 * pebblePos_[i - 1].Z);
        }
        return globalPos;
    }

    std::vector <Mdouble> getPebbleRadii(){ return pebbleRadius_; }

    // Methods setting and getting some extra Boolean properties

    // check if particle is in "Dzhanibekov" state
    bool getDzhanibekovParticle()
    {
    	return DzhanibekovParticle_;
    }

    // check if particle is "vertically oriented"
    bool getVerticallyOriented()
    {
    	return VerticallyOriented_;
    }

    // set the "Dzhanibekov" state
    void setDzhanibekovParticle( bool d)
    {
    	DzhanibekovParticle_ = d;
    }

    // set "vertically oriented" state
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

    int nPebble_;               // Number of pebbles
    
    bool DzhanibekovParticle_;  // This property is needed to quantify Dzhanibekov gas properties

    bool VerticallyOriented_;   // This property is useful for mechnical stability simulations (Gomboc, Dominos)

    Vec3D angularAcceleration_; // Clump angular acceleration

    Mdouble clumpMass_;         // Clump mass

    Mdouble viscousDamping_;    // Viscous damping parameter - for extra clump damping

    MatrixSymmetric3D clumpInertia_;        // Clump tensor of inertia

    MatrixSymmetric3D clumpInitInertia_;    // Clump initial tensor of inertia

    std::vector<Vec3D> pebblePos_;          // positions of pebbles in a local coordinate system

    std::vector<Mdouble> pebbleRadius_;     // radii of pebbles

    Matrix3D principalDirections_;          // clump's principal directions

    Matrix3D initPrincipalDirections_;      // clump's initial principal directions

    std::vector<ClumpParticle*> pebbleParticles_; // pointers to pebbles

    //Helper functions

    // Converts a Matrix3D into a MatrixSymmetric3D
    MatrixSymmetric3D MtoS( Matrix3D M){ return MatrixSymmetric3D(M.XX, M.XY, M.XZ, M.YY, M.YZ, M.ZZ);}

    // Converts a MatrixSymmetric3D into a Matrix3D
    Matrix3D StoM( MatrixSymmetric3D M){ return Matrix3D(M.XX, M.XY, M.XZ, M.XY, M.YY, M.YZ, M.XZ, M.YZ, M.ZZ);}

    // Transposes the matrix
    Matrix3D transpose(Matrix3D M){ return Matrix3D(M.XX, M.YX, M.ZX, M.XY, M.YY, M.ZY, M.XZ, M.YZ, M.ZZ);}
};

#endif
