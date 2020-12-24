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

#include "ScrewsymmetricIntersectionOfWalls.h"
#include "Particles/BaseParticle.h"
#include "WallHandler.h"
#include "DPMBase.h"
using mathsFunc::square;

ScrewsymmetricIntersectionOfWalls::ScrewsymmetricIntersectionOfWalls()
{
    logger(ERROR, "ScrewsymmetricIntersectionOfWalls() is not yet completely implemented");
    logger(DEBUG, "ScrewsymmetricIntersectionOfWalls() finished");
}

/*!
 * \param[in] other The ScrewsymmetricIntersectionOfWalls that must be copied.
 */
ScrewsymmetricIntersectionOfWalls::ScrewsymmetricIntersectionOfWalls(const ScrewsymmetricIntersectionOfWalls& other)
        : IntersectionOfWalls(other)
{
    logger(ERROR, "ScrewsymmetricIntersectionOfWalls() is not yet completely implemented");
    logger(DEBUG, "ScrewsymmetricIntersectionOfWalls(const ScrewsymmetricIntersectionOfWalls &p) finished");
}

ScrewsymmetricIntersectionOfWalls::ScrewsymmetricIntersectionOfWalls(Vec3D position, Vec3D orientation,
                                                                 std::vector<ScrewsymmetricIntersectionOfWalls::normalAndPosition> walls,
                                                                 const ParticleSpecies* species)
        : IntersectionOfWalls(walls, species)
{
    logger(ERROR, "ScrewsymmetricIntersectionOfWalls() is not yet completely implemented");
    setPosition(position);
    setOrientationViaNormal(orientation);
}

ScrewsymmetricIntersectionOfWalls::~ScrewsymmetricIntersectionOfWalls()
{
    logger(DEBUG, "~ScrewsymmetricIntersectionOfWalls() finished.");
}

/*!
 * \param[in] other The ScrewsymmetricIntersectionOfWalls that must be copied.
 */
ScrewsymmetricIntersectionOfWalls&
ScrewsymmetricIntersectionOfWalls::operator=(const ScrewsymmetricIntersectionOfWalls& other)
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
ScrewsymmetricIntersectionOfWalls* ScrewsymmetricIntersectionOfWalls::copy() const
{
    return new ScrewsymmetricIntersectionOfWalls(*this);
}

/*!
 * \details First, the particle is translated by the vector position_, then the 
 * distance normal and tangential to the orientation is computed. This normal 
 * and tangential direction is interpreted as the x and z coordinate. With the 
 * particle shifted into the XZ plane, the distance and normal is computed, as 
 * if the ScrewsymmetricIntersectionOfWalls would be a simple IntersectionOfWalls.
 * Finally, the object and the normal is rotated back to the original position.
 * 
 * See also ScrewsymmetricIntersectionOfWalls for details.
 */
bool ScrewsymmetricIntersectionOfWalls::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // squared radial position of the particle
    const Vec3D positionLabFrame = p.getPosition() - getPosition();
    const Mdouble rho2 = positionLabFrame.X * positionLabFrame.X + positionLabFrame.Y * positionLabFrame.Y;

    // should thickness be added here?
    const Mdouble wallInteractionRadius = p.getWallInteractionRadius(this) + thickness_;

    // if the particle is outside the cylinder containing the screw there is no collision
    if (//rho2 > square(rMax_ + wallInteractionRadius) ||
        positionLabFrame.Z > length_ + p.getWallInteractionRadius(this) ||
        positionLabFrame.Z < -p.getWallInteractionRadius(this))
    {
        return false;
    }
    Vec3D normalVector;
    Vec3D radialVector;
    Mdouble deltaN;
    // positionLabFrame = normalVector + radialVector + axialVector
    computeNormalRadialDeltaN(positionLabFrame, normalVector, radialVector, deltaN);

    // radial position of the particle
    const Mdouble r = std::sqrt(rho2);

    //determine wall distance, normal and contact in axissymmetric coordinates
    //and transform from axisymmetric coordinates
    Vec3D normalRN;
    if (!IntersectionOfWalls::getDistanceAndNormal(Vec3D(r, 0, deltaN), p.getWallInteractionRadius(this), distance, normalRN))
    {
        //if not in contact
        return false;
    }
    else
    {
        //if in contact
        if (r != 0)
            radialVector /= r;
        else //in this case the tangential vector is irrelevant
            logger(WARN, "Warning: Particle % is exactly on the symmetry axis of wall %", p.getIndex(), getIndex());
        //normal_Return = normalVector.Z * deltaN + normalVector.X * r;
        normal_return = normalRN.Z * normalVector + normalRN.X * radialVector;
        normal_return.normalise();
        return true;
    }
}

/*!
 * Computes the normal and radial vector for the screw at the position of the particle. Furthermore, it also computes
 *  deltaN, which is component of the particle-blade centre distance normal to the blade surface
 * \param[in] p particle whose position is used to compute the vectors and distance to the blade
 * \param[out] normalVector the vector in (x,y,z)-direction that points perpendicular to the screw blade
 * \param[out] radialVector the vector in (x,y)-direction that points outward from shaft to particle
 * \param[out] deltaN component of the particle-blade_center distance normal to the blade surface
 */
void ScrewsymmetricIntersectionOfWalls::computeNormalRadialDeltaN(const Vec3D& positionLabFrame, Vec3D& normalVector,
                                                             Vec3D& radialVector, Mdouble& deltaN) const
{
    // radial position of the particle
    const Mdouble rho = std::sqrt(positionLabFrame.X * positionLabFrame.X
                                  + positionLabFrame.Y * positionLabFrame.Y);

    // The rescaled length of the screw (length_/(2*pi*numberOfTurns_)).
    const Mdouble h = 0.5 * pitch_ / constants::pi;

    // The pitch length of the screw (length_/numberOfTurns_).
    const Mdouble deltaZ = computeDeltaZ(positionLabFrame, h, pitch_);


    // trigonometric functions relative to the particle angle
    const Mdouble cosXi = positionLabFrame.X / rho;
    const Mdouble sinXi = positionLabFrame.Y / rho;
    if (rightHandedness_)
    {
        normalVector.X = h * sinXi;
        normalVector.Y = -h * cosXi;
        normalVector.Z = rho;
    }
    else
    {
        normalVector.X = -h * cosXi;
        normalVector.Y = h * sinXi;
        normalVector.Z = rho;
    }

    // takes the right direction (+/-) of the vector and normalizes
    normalVector *= -deltaZ;
    normalVector /= normalVector.getLength();

    // radial vector at the particle position
    radialVector.X = cosXi;
    radialVector.Y = sinXi;
    radialVector.Z = 0.0;

    // The half-thickness of the screw.
    const Mdouble delta = 0.5 * thickness_;

    // cosine of the helix angle at the particle position
    const Mdouble cosEta = 1.0 / std::sqrt(1.0 + (h * h / rho / rho));

    // component of the particle-blade_center distance normal to the blade surface
    deltaN = fabs(deltaZ) * cosEta - delta;
}

/*!
 * Auxiliary function that computes the oriented axial distance between the particle's centre and the blade centre
 * \param[in] positionLabFrame distance between the starting point of the screw, start_, and the particle position
 * (Vec3D)
 * \param[in] h The rescaled length of the screw (length_/(2*pi*numberOfTurns_)) (Mdouble > 0)
 * \param[in] pitch The pitch length of the screw (length_/numberOfTurns_) (Mdouble > 0)
 * @return oriented axial distance between the particle's centre and the blade centre
 */
Mdouble ScrewsymmetricIntersectionOfWalls::computeDeltaZ(const Vec3D& positionLabFrame, Mdouble h, Mdouble pitch) const
{
    // oriented axial distance between the particle's centre and the blade centre
    Mdouble deltaZ;

    // normal to the blade at the particle position
    if (rightHandedness_) // right-handed thread
    {
        // angular coordinate of the particle
        // IMPORTANT: this angle needs to be defined in the interval [0, +2*pi[ radians!
        Mdouble xi = atan2(positionLabFrame.Y, positionLabFrame.X);
        if (xi < 0.0)
            xi += 2.0 * constants::pi;

        deltaZ = fmod(positionLabFrame.Z - h * (xi + getOrientation().getAxis().Z) -
                      static_cast<int> (positionLabFrame.Z / pitch), pitch);
        logger(DEBUG, "xi: %, deltaZ: %", xi, deltaZ);
    }
    else // left-handed thread
    {
        // angular coordinate of the particle
        // IMPORTANT: this angle needs to be defined in the interval [0, +2*pi[ radians!
        Mdouble xi = atan2(-positionLabFrame.Y, positionLabFrame.X);
        if (xi < 0.0)
            xi += 2.0 * constants::pi;
        xi += 0.5 * constants::pi;
        xi = fmod(xi, 2.0 * constants::pi);

        deltaZ = fmod(positionLabFrame.Z - h * (xi + 0.5 * constants::pi - getOrientation().getAxis().Z) -
                      static_cast<int> (positionLabFrame.Z / pitch), pitch);
        logger(DEBUG, "xi: %, deltaZ: %", xi, deltaZ);
    }

    if (deltaZ > 0.5 * pitch)
    {
        deltaZ -= pitch;
    }
    else if (deltaZ < -0.5 * pitch)
    {
        deltaZ += pitch;
    }
    return deltaZ;
}

/*!
 * \param[in] is The input stream from which the ScrewsymmetricIntersectionOfWalls
 * is read, usually a restart file.
 */
void ScrewsymmetricIntersectionOfWalls::read(std::istream& is)
{
    IntersectionOfWalls::read(is);
}

/*!
 * \param[in] os The output stream where the ScrewsymmetricIntersectionOfWalls must be written
 *  to, usually a restart file.
 */
void ScrewsymmetricIntersectionOfWalls::write(std::ostream& os) const
{
    IntersectionOfWalls::write(os);
}

/*!
 * \return The string "ScrewsymmetricIntersectionOfWalls".
 */
std::string ScrewsymmetricIntersectionOfWalls::getName() const
{
    return "ScrewsymmetricIntersectionOfWalls";
}

void ScrewsymmetricIntersectionOfWalls::setAxis(Vec3D a)
{
    setOrientationViaNormal(a);
}

void ScrewsymmetricIntersectionOfWalls::convertLimits(Vec3D& min, Vec3D& max) const
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

void ScrewsymmetricIntersectionOfWalls::writeVTK(VTKContainer& vtk) const
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
