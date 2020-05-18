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

#include <limits>

#include "RestrictedWall.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"

RestrictedWall::RestrictedWall()
{
    wall_ = nullptr;
    restriction_ = nullptr;
    logger(DEBUG, "RestrictedWall::RestrictedWall ) finished");
}

/*!
 * \param[in] w RestrictedWall that has to be copied.
 * \details First copy the attributes of the BaseWall, then copy the ones that are
 * specific for the RestrictedWall.
 */
RestrictedWall::RestrictedWall(const RestrictedWall& w)
        : BaseWall(w)
{
    wall_ = nullptr;
    restriction_ = nullptr;
    set(w.wall_, w.restriction_);
    logger(DEBUG, "RestrictedWall::RestrictedWall(const RestrictedWall &p) finished");
}

RestrictedWall::RestrictedWall(BaseWall* wall, InfiniteWall* restriction)
{
    wall_ = nullptr;
    restriction_ = nullptr;
    set(wall, restriction);
}

RestrictedWall::~RestrictedWall()
{
    if (wall_ != nullptr)
    {
        delete wall_;
        wall_ = nullptr;
    }
    if (restriction_ != nullptr)
    {
        delete restriction_;
        restriction_ = nullptr;
    }
    logger(DEBUG, "RestrictedWall::~RestrictedWall finished");
}

/*!
 * Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
 */
RestrictedWall* RestrictedWall::copy() const
{
    return new RestrictedWall(*this);
}

/*
 * \param[in] normal A Vec3D that represents the normal to the wall.
 * \param[in] point A Vec3D which is a point on the wall.
 * \details Sets the wall such that for all points x on the wall it holds that 
 * normal*x=normal*point.
 */
///\todo TW maybe the Restricted wall should be templated with the wall type such that we don't need to use new and delete.
void RestrictedWall::set(BaseWall* wall, InfiniteWall* restriction)
{
    if (wall_ != nullptr)
    {
        delete wall_;
        wall_ = nullptr;
    }
    wall_ = wall->copy();
    setSpecies(wall_->getSpecies());
    
    if (restriction_ != nullptr)
    {
        delete restriction_;
        restriction_ = nullptr;
    }
    restriction_ = restriction->copy();
    
    // std::cout << *this << std::endl;
}

/*!
 * \param[in] p BaseParticle for which the distance to the wall must be computed.
 * \param[out] distance Distance between the particle and the wall.
 * \param[out] normal_return The normal of this wall, will only be set if there is a collision.
 * \return A boolean value for whether or not there is a collision.
 * \details First the distance is checked. If there is no collision, this
 * function will return false and give the distance. If there is a collision, the
 * function will return true and give the distance and the normal vector of this wall.
 * Since this function should be called before calculating any 
 * Particle-Wall interactions, it can also be used to set the normal vector in 
 * case of curved walls.
 */
bool RestrictedWall::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    if (restriction_->getDistance(p.getPosition()) < p.getWallInteractionRadius(this))
        return wall_->getDistanceAndNormal(p, distance, normal_return);
    else
        return false;
}

/*!
 * \param[in] is The input stream from which the RestrictedWall is read.
 */
void RestrictedWall::read(std::istream& is)
{
    wall_->read(is);
    restriction_->read(is);
}

/*!
 * \param[in] os The output stream the RestrictedWall is written to.
 */
void RestrictedWall::write(std::ostream& os) const
{
    //todo check
    os << getName() << ' ';
    wall_->write(os);
    os << ' ';
    restriction_->write(os);
}

/*!
 * \return The string "RestrictedWall", which is the name of this class.
 */
std::string RestrictedWall::getName() const
{
    return "RestrictedWall";
}

/*!
 * \param[in] p Pointer to the BaseParticle which we want to check the interaction for.
 * \param[in] timeStamp The time at which we want to look at the interaction.
 * \param[in] interactionHandler A pointer to the InteractionHandler in which the interaction can be found.
 * \return A pointer to the BaseInteraction that happened between this RestrictedWall
 * and the BaseParticle at the timeStamp.
 */
///\todo Shouldn't this function be defined in BaseWall?
BaseInteraction*
RestrictedWall::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
{
    if (restriction_->getDistance(p->getPosition()) < p->getWallInteractionRadius(this))
    {
        ///\todo{setting the index of the wall is necessary to get the right index reported in fstat; however, the better way would be to make setIndex virtual.}
        wall_->setIndex(getIndex());
        return wall_->getInteractionWith(p, timeStamp, interactionHandler);
    }
    else
    {
        return nullptr;
    }
}

void RestrictedWall::writeVTK(VTKContainer& vtk) const
{
    const size_t size0 = vtk.triangleStrips.size();
    //write full BaseWall
    wall_->writeVTK(vtk);
    
    //copy out the new triangle strips, so we can modify them
    std::vector<std::vector<double>> myTriangleStrips;
    while (vtk.triangleStrips.size() > size0)
    {
        myTriangleStrips.push_back(vtk.triangleStrips.back());
        vtk.triangleStrips.pop_back();
    }
    
    //now remove points that are in restricted area
    for (const std::vector<double>& myTriangleStrip: myTriangleStrips)
    {
        unsigned counter = 0;
        std::vector<double> cell;
        //for each triangle
        for (const double c: myTriangleStrip)
        {
            Mdouble distance = restriction_->getDistance(vtk.points[c]);
            if (distance < 0)
            {
                //if the current point is inside the restricted volume, write it
                cell.push_back(c);
                counter = 0;
            }
            else
            {
                if (counter >= 2)
                {
                    //if the current point is not in restriction, don't write it and flush the cell into the cells array
                    if (cell.size() > 2)
                        vtk.triangleStrips.push_back(cell);
                    cell.clear();
                    counter = 0;
                }
                vtk.points[c] += distance * restriction_->getNormal();
                cell.push_back(c);
                counter++;
            }
        }
        if (cell.size() > 2)
        {
            //if the current point is not in restriction, don't write it and flush the cell into the cells array
            vtk.triangleStrips.push_back(cell);
            cell.clear();
        }
    }
}
