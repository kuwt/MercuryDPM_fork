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

#ifndef BASEWALL_H
#define BASEWALL_H

#include "BaseInteractable.h"
#include "Particles/SuperQuadricParticle.h"

class WallHandler;

class BaseParticle;

struct VTKContainer
{
    std::vector<Vec3D> points;
    std::vector<std::vector<double>> triangleStrips;
};

/*!
 * \brief Basic class for walls.
 * \details Class from which all walls inherit. Please note the getVelocity can 
 * for some walls be dependent on which point on the wall is meant.
 */
class BaseWall : public BaseInteractable
{
public:
    /*!
     * \brief Default constructor.
     */
    BaseWall();
    
    /*!
     * \brief Copy constructor.
     */
    BaseWall(const BaseWall& w);
    
    /*!
     * \brief Default destructor. 
     */
    ~BaseWall() override;
    
    /*!
     * \brief Pure virtual function that can be overwritten in inherited classes in order to copy a BaseWall.
     * \return A pointer to the new BaseWall.
     */
    virtual BaseWall* copy() const = 0;
    
    /*!
     * \brief Function that reads a BaseWall from an input stream, usually a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Function that writes a BaseWall to an output stream, usually a restart file.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Pure virtual function that computes the distance of a BaseParticle to this wall and returns the normal of this wall if there is a collision.
     * \details Beware, the distance and normal are output parameters, not return values!
     * \param[in] P Reference to the BaseParticle we want to compute the distance to the BaseWall of.
     * \param[out] distance Distance of the BaseParticle to the BaseWall.
     * \param[out] normal_return The normal of the wall. Is only given if there is a collision.
     * \return A boolean which indicates if there is a collision between the BaseParticle and the wall.
     */
    virtual bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const = 0;
    
    
    virtual bool
    getDistanceNormalOverlap(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return, Mdouble& overlap) const;
    
    virtual bool getDistanceNormalOverlapSuperquadric(const SuperQuadricParticle& p, Mdouble& distance, Vec3D& normal_return,
                                                      Mdouble& overlap) const;
    
    virtual Vec3D
    getFurthestPointSuperQuadric(const Vec3D& normalBodyFixed, const Vec3D& axes, Mdouble eps1, Mdouble eps2) const;
    
    /*!
     * \brief A function which sets the WallHandler for this BaseWall.
     */
    virtual void setHandler(WallHandler* handler);
    
    /*!
     * \brief A function which returns the WallHandler that handles this BaseWall.
     */
    WallHandler* getHandler() const;
    
    /*!
     * \deprecated TW: this function should be taken out and replaced by setSpecies
     * \brief Define the species of this wall using the index of the species in the SpeciesHandler in this DPMBase.
     */
    void setIndSpecies(unsigned int indSpecies) override;
    
    /*!
     * \brief Defines the species of the current wall.
     *
     */
    void setSpecies(const ParticleSpecies* species);
    
    /*!
     * \details Sets all walls (unlike particles) to be inherently fixed - i.e. the bool "is fixed" will by definition
     * return "true" when called for a wall (i.e. any BaseWall tyope object).
     */
    bool isFixed() const override;

    /*! 
     * \brief Slowly adjusts the force on a wall towards a specified goal, by adjusting (prescribing) the velocity of the wall.
     */
    void setForceControl(Vec3D forceGoal, Vec3D gainFactor, Vec3D baseVelocity={0,0,0});
    
    /*!
     * if isLocal returns true and the DPM class derived from MercuryBase, the hGrid will be used to find wall-particle contacts, using min/max.
     */
    virtual bool isLocal(Vec3D& min, Vec3D& max) const
    { return false; }
    
    // returns the point intersecting a wall (x-p).n=0 and a line x=p0+t(p1-p0)
    bool getLinePlaneIntersect(Vec3D& intersect, const Vec3D& p0, const Vec3D& p1, const Vec3D& n, const Vec3D& p);
    
    /*!
     * Checks if point is in wall (if close to the wall, the point is assumed out of the wall)
    */bool isInsideWallVTK(const Vec3D& point, const Vec3D& normal, const Vec3D& position) const;
    
    /*!
     * Moves point0 onto the surface of the wallsuch that the direction of the edge from point0 to point1 is preserved.
    */
    void projectOntoWallVTK(Vec3D& point0, const Vec3D& point1, const Vec3D& normal, const Vec3D& position) const;
    
    /*!
     * Intersects a of set VTK points representing a wall with a half-space defined by position and normal;
     * This is used in createVTK to restrict the VTK points representing a wall to the inside of the domain.
*/
    void intersectVTK(std::vector<Vec3D>& points, Vec3D normal, Vec3D position) const;
    
    
    virtual BaseInteraction*
    getInteractionWithSuperQuad(SuperQuadricParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler);
    
    /*!
     * adds extra information to the points and triangleStrips vectors needed to plot the wall in vtk format
     * @param points Coordinates of the vertices of the triangulated surfaces (in the VTK file this is called POINTS)
     * @param triangleStrips Indices of three vertices forming one triangulated surface (in the VTK file this is called CELL)
     */
    virtual void writeVTK(VTKContainer& vtk) const;
    
    void getVTK(std::vector<Vec3D>& points, std::vector<std::vector<double>>& triangleStrips)
    {}
    
    const Vec3D getAxis() const;
    
    ///\brief Returns the interaction between this wall and a given particle, nullptr if there is no interaction.
    BaseInteraction*
    getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) override;
    
    bool getVTKVisibility() const;
    
    void setVTKVisibility(bool vtkVisibility);
    
    /*!
     * \brief Takes the points provided and adds a triangle strip connecting these points to the vtk container
     */
    static void addToVTK(const std::vector<Vec3D>& points, VTKContainer& vtk);
    
    //todo how do i write a copy and add function?
    void addRenderedWall(BaseWall* w);

    BaseWall* getRenderedWall(size_t i) const;

    void removeRenderedWalls();

    void renderWall(VTKContainer& vtk);

    void addParticlesAtWall(unsigned numElements = 50);

    void setVelocityControl(Vec3D forceGoal, Vec3D gainFactor, Vec3D baseVelocity);

private:
    /*!
     * A pointer to the WallHandler that handles this BaseWall.
     */
    WallHandler* handler_; ///
    
    bool vtkVisibility_ = true;
    
    /*!
     * A vector of walls that gets rendered instead of the actual wall
     */
    std::vector<BaseWall*> renderedWalls_;

};

#endif
