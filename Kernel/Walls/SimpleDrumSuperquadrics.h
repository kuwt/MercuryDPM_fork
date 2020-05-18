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

#ifndef AXISYMMETRICWALL_H
#define AXISYMMETRICWALL_H

#include "Math/Vector.h"
#include "Walls/InfiniteWall.h"

/*!
 * \brief A drum in xz-direction with centre at the origin with a certain radius. Usable with superquadric particles
 * \details Notice that this wall contains a mutable infinite wall. This is because for contact-detection we map the
 * infinite wall on a point with the correct distance from the origin, thereby linearising the drum. Since this infinite
 * wall has to change normal and position for every interaction, it must be changeable in the const method
 * getDistanceAndNormal and therefore the infinite wall must be mutable (or volatile, but that is not necessary here).
 *
 * \todo Document
 * \todo Definitions in source-file
 * \todo VTK-write function
 * \todo Test
 */
class SimpleDrumSuperquadrics : public BaseWall
{
public:
    /*!
     * \brief Default constructor.
     */
    SimpleDrumSuperquadrics();
    
    /*!
     * \brief Copy constructor.
     */
    SimpleDrumSuperquadrics(const SimpleDrumSuperquadrics& p);
    
    
    /*!
     * \brief Destructor.
     */
    ~SimpleDrumSuperquadrics() override;
    
    /*!
     * \brief Copy assignment operator.
     */
    SimpleDrumSuperquadrics& operator=(const SimpleDrumSuperquadrics& other);
    
    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    SimpleDrumSuperquadrics* copy() const final;
    
    /*!
     * \brief Computes the distance from the wall for a given BaseParticle and 
     * returns true if there is a collision. If there is a collision, also 
     * return the normal vector.
     */
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const override;
    
    bool getDistanceNormalOverlapSuperquadric(const SuperQuadricParticle& p, Mdouble& distance, Vec3D& normal_return,
                                              Mdouble& overlap) const override;
    
    /*!
     * \brief reads wall
     */
    void read(std::istream& is) final;
    
    /*!
     * \brief outputs wall
     */
    void write(std::ostream& os) const final;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const final;
    
    void setAxis(Vec3D a);
    
    void setRadius(const Mdouble& radius)
    {
        radius_ = radius;
        logger(INFO, "setting drum radius to %", radius_);
    }
    
    BaseInteraction* getInteractionWithSuperQuad(SuperQuadricParticle* p, unsigned timeStamp,
                                                                        InteractionHandler* interactionHandler) override
    {
        return wall.getInteractionWithSuperQuad(p, timeStamp, interactionHandler);
    }
    
    Vec3D getFurthestPointSuperQuadric(const Vec3D& normalBodyFixed, const Vec3D& axes, Mdouble eps1, Mdouble eps2) const override
    {
        return wall.getFurthestPointSuperQuadric(normalBodyFixed, axes, eps1, eps2);
    }

private:
    mutable InfiniteWall wall;
    Mdouble radius_;
};


#endif
