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

#include "HorizontalBaseScrew.h"
#include "Particles/BaseParticle.h"
#include "WallHandler.h"
#include "DPMBase.h"

HorizontalBaseScrew::HorizontalBaseScrew()
{
    //BEGIN CHANGE
    sinA2Max_ = 0;
    timeMin_ = 0;
    //END CHANGE
    logger(DEBUG, "HorizontalBaseScrew() finished");
}

/*!
 * \param[in] other The HorizontalBaseScrew that must be copied.
 */
HorizontalBaseScrew::HorizontalBaseScrew(const HorizontalBaseScrew &other)
        : IntersectionOfWalls(other)
{
    //BEGIN CHANGE
    sinA2Max_ = other.sinA2Max_;
    timeMin_ = other.timeMin_;
    //END CHANGE
    logger(DEBUG, "HorizontalBaseScrew(const HorizontalBaseScrew &p) finished");
}

HorizontalBaseScrew::HorizontalBaseScrew(Vec3D position, Vec3D orientation, std::vector<HorizontalBaseScrew::normalAndPosition> walls, const ParticleSpecies* species, Mdouble sinA2Max, Mdouble timeMin)
        : IntersectionOfWalls(walls, species), sinA2Max_(sinA2Max), timeMin_(timeMin)
{
    setPosition(position);
    setOrientationViaNormal(orientation);
}

HorizontalBaseScrew::~HorizontalBaseScrew()
{
    logger(DEBUG, "~HorizontalBaseScrew() finished.");
}

/*!
 * \param[in] other The HorizontalBaseScrew that must be copied.
 */
HorizontalBaseScrew& HorizontalBaseScrew::operator=(const HorizontalBaseScrew& other)
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
HorizontalBaseScrew* HorizontalBaseScrew::copy() const
{
    return new HorizontalBaseScrew(*this);
}

/*!
 * \details First, the particle is translated by the vector position_, then the 
 * distance normal and tangential to the orientation is computed. This normal 
 * and tangential direction is interpreted as the x and z coordinate. With the 
 * particle shifted into the XZ plane, the distance and normal is computed, as 
 * if the HorizontalBaseScrew would be a simple IntersectionOfWalls.
 * Finally, the object and the normal is rotated back to the original position.
 * 
 * See also HorizontalBaseScrew for details.
 */
bool HorizontalBaseScrew::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normalReturn) const
{
    //transform to axisymmetric coordinates
    //move the coordinate system to the axis origin, so pOrigin=(xhat,yhat,zhat)
    Vec3D pOrigin = p.getPosition() -getPosition();
    ///\todo check, maybe orientation has to be normalised differently for axisymmetric walls (or axis needs to be normalised)
    Vec3D axis = getOrientation().getAxis();
    Mdouble a = Vec3D::dot(pOrigin, axis);
    //Vec3D(r,a) is the projection into the radial-axial plane
    Vec3D radialDirection = pOrigin - a * axis;
    Mdouble r = radialDirection.getLength();

    //BEGIN CHANGE
    //if in contact
    if (r!=0)
        radialDirection /= r;
    else //in this case the tangential vector is irrelevant
        logger(WARN, "Warning: Particle % is exactly on the symmetry axis of wall %", p.getIndex(), getIndex());

    //tangential is the projection into the (xhat,yhat) plane
    Mdouble angle = getAngularVelocity().Z*(getHandler()->getDPMBase()->getTime()-timeMin_);
    Vec3D orientation = Vec3D(cos(angle),sin(angle),0);
    Mdouble sinA2 = Vec3D::getLengthSquared(Vec3D::cross(radialDirection,orientation));
    if (sinA2>sinA2Max_) {
        a+=sinA2-sinA2Max_;
    }
    //END CHANGE

    //determine wall distance, normal and contact in axissymmetric coordinates
    //and transform from axisymmetric coordinates
    Vec3D normal;
    if (!IntersectionOfWalls::getDistanceAndNormal(Vec3D(r, 0, a), p.getWallInteractionRadius(this), distance, normal))
    {
        //if not in contact
        return false;
    }
    else
    {
        //if in contact
        if (r!=0)
            radialDirection /= r;
        else //in this case the tangential vector is irrelevant
            logger(WARN, "Warning: Particle % is exactly on the symmetry axis of wall %", p.getIndex(), getIndex());
        normalReturn = normal.Z * axis + normal.X * radialDirection;
        //logger.assert(normalReturn.Y==0,"Error");
        return true;
    }
}

/*!
 * \param[in] is The input stream from which the HorizontalBaseScrew
 * is read, usually a restart file.
 */
void HorizontalBaseScrew::read(std::istream& is)
{
    IntersectionOfWalls::read(is);
}

/*!
 * \param[in] os The output stream where the HorizontalBaseScrew must be written
 *  to, usually a restart file.
 */
void HorizontalBaseScrew::write(std::ostream& os) const
{
    IntersectionOfWalls::write(os);
}

/*!
 * \return The string "HorizontalBaseScrew".
 */
std::string HorizontalBaseScrew::getName() const
{
    return "HorizontalBaseScrew";
}

void HorizontalBaseScrew::setAxis(Vec3D a)
{
    setOrientationViaNormal(a);
}

void HorizontalBaseScrew::convertLimits (Vec3D& min, Vec3D& max) const {
    //define the box in the rz plane 
    double x[2] = {min.X, max.X};
    double y[2] = {min.Y, max.Y};
    double z[2] = {min.Z, max.Z};
    Quaternion q = getOrientation();
    Vec3D o = q.getAxis();

    //look at the corner points of the domain and see that the wall is displayed in the domain (possibly a bit outside as well)
    min -= getPosition();
    q.rotate(min); //set min/max initial values to values of first corner point
    max = min;
    Mdouble maxXSquared = 0;
    //
    for (unsigned i=0; i<2; i++)
        for (unsigned j=0; j<2; j++)
            for (unsigned k=0; k<2; k++) {
                Vec3D p = Vec3D(x[i],y[j],z[k])-getPosition();
                Mdouble Z = Vec3D::dot(p,o);
                Mdouble XSquared = Vec3D::getLengthSquared(p-Z*o);
                if (XSquared>maxXSquared) maxXSquared=XSquared;
                if (Z<min.Z) min.Z=Z;
                if (Z>max.Z) max.Z=Z;
            }
    min.X = 0;
    max.X = sqrt(maxXSquared);
    min.Y = 0;
    max.Y = 0.001; //unused
    //logger(ERROR,"min % max %",min,max);
}

void HorizontalBaseScrew::writeVTK (VTKContainer& vtk) const
{
    for (auto wall=wallObjects_.begin(); wall!=wallObjects_.end(); wall++)
    {
        //plot each of the intersecting walls
        std::vector<Vec3D> myPoints;

        //first create a slice of non-rotated wall in the xz plane, 0<y<1            
        Vec3D min = getHandler()->getDPMBase()->getMin();
        Vec3D max = getHandler()->getDPMBase()->getMax();
        convertLimits(min, max);

        //create the basic slice for the first wall using the InfiniteWall routine
        wall->createVTK (myPoints,min,max);

        //create intersections with the other walls, similar to the IntersectionOfWalls routine
        for (auto other=wallObjects_.begin(); other!=wallObjects_.end(); other++)
        {
            if (other!=wall)
            {
                intersectVTK(myPoints, -other->getNormal(), other->getPosition());
            }
        }

        //only keep the y=0 values
        std::vector<Vec3D> rzVec;
        for (auto& p: myPoints)
        {
            if (p.Y==0)
            {
                rzVec.push_back(p);
            }
        }
        if (rzVec.size()==0)
            return;

        //create points on the unit circle
        unsigned nr = 180;
        //BEGIN CHANGE
        std::vector<Vec3D> xyVec;
        Mdouble a1 = getAngularVelocity().Z*(getHandler()->getDPMBase()->getTime()-timeMin_);
        for (unsigned ir=0; ir<nr; ir++) {
            Mdouble angle = 2.0 * constants::pi * ir/nr;

            Mdouble sinA2 = sin(a1 - angle); sinA2 *= sinA2;
            Mdouble dSin = sinA2-sinA2Max_;

            xyVec.push_back(Vec3D(cos(angle),sin(angle),fmax(dSin,0)));
        }
        //END CHANGE

        //now create rings of points on the axisym. shape
        ///\bug once the quaternions are implemented, we can orient these walls properly
        unsigned nPoints = vtk.points.size();
        Vec3D p;
        Vec3D o = getOrientation().getAxis();
        for (auto rz : rzVec)
        {
            for (auto xy : xyVec)
            {
                //BEGIN CHANGE
                p = Vec3D(rz.X * xy.X, rz.X * xy.Y,fmax(rz.Z - xy.Z,0.0));
                //getOrientation().rotate(p);
                //END CHANGE
                p += getPosition();
                vtk.points.push_back(p);
            }
        }

        //finally create the connectivity matri to plot shell-like triangle strips.
        unsigned nz = rzVec.size();
        unsigned nCells = vtk.triangleStrips.size();
        for (unsigned iz=0; iz<nz-1; iz++) {
            std::vector<double> cell;
            cell.reserve(2*nr+2);
            for (unsigned ir=0; ir<nr; ir++) {
                cell.push_back(nPoints+ir+iz*nr);
                cell.push_back(nPoints+ir+(iz+1)*nr);
            }
            cell.push_back(nPoints+iz*nr);
            cell.push_back(nPoints+(iz+1)*nr);
            vtk.triangleStrips.push_back(cell);
        }
    }

}
