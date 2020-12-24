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

#include "AxisymmetricIntersectionOfWalls.h"
#include "Particles/BaseParticle.h"
#include "WallHandler.h"
#include "DPMBase.h"

AxisymmetricIntersectionOfWalls::AxisymmetricIntersectionOfWalls()
{
    logger(DEBUG, "AxisymmetricIntersectionOfWalls() finished");
}

/*!
 * \param[in] other The AxisymmetricIntersectionOfWalls that must be copied.
 */
AxisymmetricIntersectionOfWalls::AxisymmetricIntersectionOfWalls(const AxisymmetricIntersectionOfWalls& other)
        : IntersectionOfWalls(other)
{
    logger(DEBUG, "AxisymmetricIntersectionOfWalls(const AxisymmetricIntersectionOfWalls &p) finished");
}

AxisymmetricIntersectionOfWalls::AxisymmetricIntersectionOfWalls(Vec3D position, Vec3D orientation,
                                                                 std::vector<AxisymmetricIntersectionOfWalls::normalAndPosition> walls,
                                                                 const ParticleSpecies* species)
        : IntersectionOfWalls(walls, species)
{
    setPosition(position);
    setOrientationViaNormal(orientation);
}

AxisymmetricIntersectionOfWalls::~AxisymmetricIntersectionOfWalls()
{
    logger(DEBUG, "~AxisymmetricIntersectionOfWalls() finished.");
}

/*!
 * \param[in] other The AxisymmetricIntersectionOfWalls that must be copied.
 */
AxisymmetricIntersectionOfWalls&
AxisymmetricIntersectionOfWalls::operator=(const AxisymmetricIntersectionOfWalls& other)
{
    if (this == &other)
    {
        return *this;
    }
    else
    {
        return *(other.copy());
    }
}

/*!
 * \return pointer to a IntersectionOfWalls object allocated using new.
 */
AxisymmetricIntersectionOfWalls* AxisymmetricIntersectionOfWalls::copy() const
{
    return new AxisymmetricIntersectionOfWalls(*this);
}

/*!
 * \details First, the particle is translated by the vector position_, then the 
 * distance normal and tangential to the orientation is computed. This normal 
 * and tangential direction is interpreted as the x and z coordinate. With the 
 * particle shifted into the XZ plane, the distance and normal is computed, as 
 * if the AxisymmetricIntersectionOfWalls would be a simple IntersectionOfWalls.
 * Finally, the object and the normal is rotated back to the original position.
 * 
 * See also AxisymmetricIntersectionOfWalls for details.
 */
bool AxisymmetricIntersectionOfWalls::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance,
                                                           Vec3D& normalReturn) const
{
    //transform to axisymmetric coordinates
    //move the coordinate system to the axis origin, so pOrigin=(xhat,yhat,zhat)
    Vec3D pOrigin = p.getPosition() - getPosition();
    ///\todo check, maybe orientation has to be normalized differently for axisymmetric walls (or axis needs to be normalized)
    Vec3D axis = getOrientation().getAxis();
    Mdouble a = Vec3D::dot(pOrigin, axis);
    //Vec3D(r,a) is the projection into the radial-axial plane
    Vec3D radialDirection = pOrigin - a * axis;
    Mdouble r = radialDirection.getLength();
    Vec3D normal;
    //determine wall distance, normal and contact in axissymmetric coordinates
    //and transform from axisymmetric coordinates
    if (!IntersectionOfWalls::getDistanceAndNormal(Vec3D(r, 0, a), p.getWallInteractionRadius(this), distance, normal))
    {
        //if not in contact
        return false;
    }
    else
    {
        //if in contact
        if (r != 0)
            radialDirection /= r;
        else //in this case the tangential vector is irrelevant
            logger(WARN, "Warning: Particle % is exactly on the symmetry axis of wall %", p.getIndex(), getIndex());
        normalReturn = normal.Z * axis + normal.X * radialDirection;
        //logger.assert(normalReturn.Y==0,"Error");
        return true;
    }
}

/*!
 * \param[in] is The input stream from which the AxisymmetricIntersectionOfWalls
 * is read, usually a restart file.
 */
void AxisymmetricIntersectionOfWalls::read(std::istream& is)
{
    IntersectionOfWalls::read(is);
}

/*!
 * \param[in] os The output stream where the AxisymmetricIntersectionOfWalls must be written
 *  to, usually a restart file.
 */
void AxisymmetricIntersectionOfWalls::write(std::ostream& os) const
{
    IntersectionOfWalls::write(os);
}

/*!
 * \return The string "AxisymmetricIntersectionOfWalls".
 */
std::string AxisymmetricIntersectionOfWalls::getName() const
{
    return "AxisymmetricIntersectionOfWalls";
}

void AxisymmetricIntersectionOfWalls::setAxis(Vec3D a)
{
    setOrientationViaNormal(a);
}

void AxisymmetricIntersectionOfWalls::convertLimits(Vec3D& min, Vec3D& max) const
{
    Quaternion q = getOrientation();
    Vec3D rMin = min - getPosition();
    q.rotateBack(rMin); //set min/max initial values to values of first corner point
    Vec3D rMax = max - getPosition();
    q.rotateBack(rMax); //set min/max initial values to values of first corner point
    
    Mdouble r = std::sqrt(std::max(rMax.Y * rMax.Y + rMax.Z * rMax.Z, rMin.Y * rMin.Y + rMin.Z * rMin.Z));
    max = Vec3D(r, 0.001, std::max(rMin.X,rMax.X));
    min = Vec3D(0, 0, std::min(rMin.X,rMax.X));
    //std::cout << "r=" << r << std::endl;
}

void AxisymmetricIntersectionOfWalls::writeVTK(VTKContainer& vtk) const
{
    for (auto wall = wallObjects_.begin(); wall != wallObjects_.end(); wall++)
    {
        //plot each of the intersecting walls
        std::vector<Vec3D> myPoints;
        
        //first create a slice of non-rotated wall in the xz plane, 0<y<1            
        Vec3D min = getHandler()->getDPMBase()->getMin();
        Vec3D max = getHandler()->getDPMBase()->getMax();
        convertLimits(min, max);
        
        //create the basic slice for the first wall using the InfiniteWall routine
        wall->createVTK(myPoints, min, max);
        
        //create intersections with the other walls, similar to the IntersectionOfWalls routine
        for (auto other = wallObjects_.begin(); other != wallObjects_.end(); other++)
        {
            if (other != wall)
            {
                intersectVTK(myPoints, -other->getNormal(), other->getPosition());
            }
        }
        
        //only keep the y=0 values
        std::vector<Vec3D> rzVec;
        for (auto& p: myPoints)
        {
            if (p.Y == 0)
            {
                rzVec.push_back(p);
            }
        }
        if (rzVec.empty())
            return;
        
        //create points on the unit circle
        unsigned nr = 180;
        struct XY
        {
            double X;
            double Y;
        };
        std::vector<XY> xyVec;
        for (unsigned ir = 0; ir < nr; ir++)
        {
            Mdouble angle = 2.0 * constants::pi * ir / nr;
            xyVec.push_back({cos(angle), sin(angle)});
        }
        
        //now create rings of points on the axisym. shape
        ///\bug once the quaternions are implemented, we can orient these walls properly
        unsigned long nPoints = vtk.points.size();
        Vec3D p;
        //Vec3D o = getOrientation().getAxis();
        for (auto rz : rzVec)
        {
            for (auto xy : xyVec)
            {
                p = Vec3D(rz.Z, rz.X * xy.X, rz.X * xy.Y);
                getOrientation().rotate(p);
                p += getPosition();
                vtk.points.push_back(p);
            }
        }
        
        //finally create the connectivity matri to plot shell-like triangle strips.
        unsigned long nz = rzVec.size();
        unsigned long nCells = vtk.triangleStrips.size();
        for (unsigned iz = 0; iz < nz - 1; iz++)
        {
            std::vector<double> cell;
            cell.reserve(2 * nr + 2);
            for (unsigned ir = 0; ir < nr; ir++)
            {
                cell.push_back(nPoints + ir + iz * nr);
                cell.push_back(nPoints + ir + (iz + 1) * nr);
            }
            cell.push_back(nPoints + iz * nr);
            cell.push_back(nPoints + (iz + 1) * nr);
            vtk.triangleStrips.push_back(cell);
        }
    }
    
}
