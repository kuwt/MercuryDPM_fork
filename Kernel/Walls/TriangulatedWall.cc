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

#include "TriangulatedWall.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

TriangulatedWall::TriangulatedWall()
{
    logger(DEBUG, "TriangulatedWall() constructed.");
}

/*!
 * \param[in] other The TriangulatedWall that must be copied.
 */
TriangulatedWall::TriangulatedWall(const TriangulatedWall& other)
        : BaseWall(other)
{
    face_ = other.face_;
    vertex_ = other.vertex_;
    //now reset the pointers
    for (unsigned f = 0; f < face_.size(); f++)
    {
        for (unsigned i = 0; i < 3; i++)
        {
            if (other.face_[f].neighbor[i])
            {
                face_[f].neighbor[i] = &face_[other.face_[f].neighbor[i] - &other.face_[0]];
            } //else nullptr
            face_[f].vertex[i] = &vertex_[other.face_[f].vertex[i] - &other.vertex_[0]];
        }
    }
    logger(DEBUG, "TriangulatedWall(TriangulatedWall&) constructed.");
}

TriangulatedWall::TriangulatedWall(std::string filename, const ParticleSpecies* species)
{
    setSpecies(species);
    readVTK(filename);
}

/**
 * \param[in] filename name of vtk input file, e.g. TriangulatedWallSelfTest.vtk
 */
void TriangulatedWall::readVTK(std::string filename)
{
    std::fstream file;
    file.open(filename.c_str(), std::ios::in);
    logger.assert_always(file.is_open(), "File opening failed: %", filename);
    
    std::string dummy;
    getline(file, dummy);
    getline(file, dummy);
    getline(file, dummy);
    getline(file, dummy);
    //read points
    unsigned num;
    file >> dummy >> num >> dummy;
    vertex_.reserve(num);
    Vec3D v;
    for (unsigned i = 0; i < num; i++)
    {
        file >> v.X >> v.Y >> v.Z;
        vertex_.push_back(v);
    }
    //read faces
    file >> dummy >> num >> dummy;
    face_.reserve(num);
    Face f;
    unsigned id0, id1, id2;
    for (unsigned i = 0; i < num; i++)
    {
        file >> dummy >> id0 >> id1 >> id2;
        f.vertex[0] = &vertex_[id0];
        f.vertex[1] = &vertex_[id1];
        f.vertex[2] = &vertex_[id2];
        face_.push_back(f);
    }
    file >> dummy;
    //set normals and positions
    for (auto& face: face_)
    {
        face.normal = Vec3D::getUnitVector(
                Vec3D::cross(*face.vertex[1] - *face.vertex[0], *face.vertex[2] - *face.vertex[0]));
    }
    //set neighbours
    for (auto face0 = face_.begin(); face0 + 1 != face_.end(); face0++)
    {
        for (auto face1 = face0 + 1; face1 != face_.end(); face1++)
        {
            if (face0->vertex[0] == face1->vertex[0])
            {
                if (face0->vertex[1] == face1->vertex[2])
                { //edge 0=2
                    face0->neighbor[0] = &*face1;
                    face1->neighbor[2] = &*face0;
                }
                else if (face0->vertex[2] == face1->vertex[1])
                { //edge 2=0
                    face0->neighbor[2] = &*face1;
                    face1->neighbor[0] = &*face0;
                }
            }
            else if (face0->vertex[0] == face1->vertex[1])
            {
                if (face0->vertex[1] == face1->vertex[0])
                { //edge 0=0
                    face0->neighbor[0] = &*face1;
                    face1->neighbor[0] = &*face0;
                }
                else if (face0->vertex[2] == face1->vertex[2])
                { //edge 2=1
                    face0->neighbor[2] = &*face1;
                    face1->neighbor[1] = &*face0;
                }
            }
            else if (face0->vertex[0] == face1->vertex[2])
            {
                if (face0->vertex[1] == face1->vertex[1])
                { //edge 0=1
                    face0->neighbor[0] = &*face1;
                    face1->neighbor[1] = &*face0;
                }
                else if (face0->vertex[2] == face1->vertex[0])
                { //edge 2=2
                    face0->neighbor[2] = &*face1;
                    face1->neighbor[2] = &*face0;
                }
            }
            else if (face0->vertex[1] == face1->vertex[0])
            {
                if (face0->vertex[2] == face1->vertex[2])
                { //edge 1=2
                    face0->neighbor[1] = &*face1;
                    face1->neighbor[2] = &*face0;
                }
            }
            else if (face0->vertex[1] == face1->vertex[1])
            {
                if (face0->vertex[2] == face1->vertex[0])
                { //edge 1=0
                    face0->neighbor[1] = &*face1;
                    face1->neighbor[0] = &*face0;
                }
            }
            else if (face0->vertex[1] == face1->vertex[2])
            {
                if (face0->vertex[2] == face1->vertex[1])
                { //edge 1=1
                    face0->neighbor[1] = &*face1;
                    face1->neighbor[1] = &*face0;
                }
            }
            
        }
    }
    //set edge normals (inwards facing)
    for (auto& face: face_)
    {
        face.edgeNormal[0] = Vec3D::getUnitVector(Vec3D::cross(face.normal, *face.vertex[1] - *face.vertex[0]));
        face.edgeNormal[1] = Vec3D::getUnitVector(Vec3D::cross(face.normal, *face.vertex[2] - *face.vertex[1]));
        face.edgeNormal[2] = Vec3D::getUnitVector(Vec3D::cross(face.normal, *face.vertex[0] - *face.vertex[2]));
    }
    file.close();
}

TriangulatedWall::~TriangulatedWall()
{
    logger(DEBUG, "~TriangulatedWall() has been called.");
}

/*!
 * \param[in] other The TriangulatedWall that must be copied.
 */
TriangulatedWall& TriangulatedWall::operator=(const TriangulatedWall& other)
{
    logger(DEBUG, "TriangulatedWall::operator= called.");
    if (this == &other)
    {
        return *this;
    }
    return *(other.copy());
}

/*!
 * \return pointer to a TriangulatedWall object allocated using new.
 */
TriangulatedWall* TriangulatedWall::copy() const
{
    return new TriangulatedWall(*this);
}

/*!
 * \param[in] p BaseParticle we want to calculate the distance and whether it collided of.
 * \param[out] distance The distance of the BaseParticle to this wall.
 * \param[out] normal_return If there was a collision, the normal vector to this wall will be placed here.
 * \return A boolean which says whether or not there was a collision.
 * \details This function computes whether or not there is a collision between 
 * a given BaseParticle and this TriangulatedWall. If there is a collision, this
 * function also computes the distance between the BaseParticle and TriangulatedWall
 * and the normal of the TriangulatedWall at the intersection point. It does
 * this by calling TriangulatedWall::getDistanceAndNormal(const Vec3D& , Mdouble , Mdouble&, Vec3D&) const.
 * Since this function should be called before calculating any 
 * Particle-Wall interactions, it can also be used to set the normal vector in 
 * case of curved walls.
 *
 * NOTE: THIS ONLY RETURNS ONE OF POSSIBLY MANY INTERACTIONS; it's only used for finding out if interactions exist, so should be fine
 */
bool TriangulatedWall::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    Mdouble interactionRadius = p.getWallInteractionRadius(this);
    //it's important to use a reference here, as the faces use pointers
    for (const auto& face : face_)
    {
        if (face.getDistanceAndNormal(p, distance, normal_return, interactionRadius))
            return true;
    }
    return false;
}

/*!
 * \param[in] move A reference to a Vec3D that denotes the direction and length 
 * it should be moved with.
 * \details A function that moves the TriangulatedWall in a certain direction
 * by both moving the walls and all intersections. Note that the directions of the 
 * intersections are not moved since they don't change when moving the TriangulatedWall
 * as a whole.
 * \todo We should use the position_ and orientation_ of the TriangulatedWall; 
 * that way, TriangulatedWall can be moved with the standard BaseInteractable::move function, 
 * getting rid of an anomaly in the code and removing the virtual from the move function. \author weinhartt
 */
void TriangulatedWall::move(const Vec3D& move)
{
    BaseInteractable::move(move);
    for (auto& v : vertex_)
    {
        v += move;
    }
}

/*!
 * \param[in] is The input stream from which the TriangulatedWall is read, usually a restart file.
 */
void TriangulatedWall::read(std::istream& is)
{
    ///\todo
}

/*!
 * \param[in] os The output stream where the TriangulatedWall must be written
 *  to, usually a restart file.
 */
void TriangulatedWall::write(std::ostream& os) const
{
    os << "Vertices " << vertex_.size();
    for (const auto& vertex: vertex_)
        os << "  " << vertex;
    os << std::endl;
    //os << "Faces " << face_.size() << std::endl;
    unsigned counter = 0;
    for (const auto& face: face_)
    {
        os << "Face " << counter++
           << " vertex " << face.vertex[0] - &vertex_[0]
           << " " << face.vertex[1] - &vertex_[0]
           << " " << face.vertex[2] - &vertex_[0]
           << " neighbor " << (face.neighbor[0] ? (face.neighbor[0] - &face_[0]) : -1)
           << " " << (face.neighbor[1] ? (face.neighbor[1] - &face_[0]) : -1)
           << " " << (face.neighbor[2] ? (face.neighbor[2] - &face_[0]) : -1)
           << " normal " << face.normal
           << " edgeNormal " << face.edgeNormal[0]
           << "  " << face.edgeNormal[1]
           << "  " << face.edgeNormal[2]
           << std::endl;
    }
}

/*!
 * \return The string "TriangulatedWall".
 */
std::string TriangulatedWall::getName() const
{
    return "TriangulatedWall";
}

/*!
 * \param[in] p Pointer to the BaseParticle which we want to check the interaction for.
 * \param[in] timeStamp The time at which we want to look at the interaction.
 * \param[in] interactionHandler A pointer to the InteractionHandler in which the interaction can be found.
 * \return A pointer to the BaseInteraction that happened between this InfiniteWall
 * and the BaseParticle at the timeStamp.
 */
BaseInteraction* TriangulatedWall::getInteractionWith(BaseParticle* p, unsigned timeStamp,
                                                                   InteractionHandler* interactionHandler)
{
    Mdouble distance;
    Vec3D normal;
    Mdouble interactionRadius = p->getWallInteractionRadius(this);
    for (const auto& face : face_)
    {
        if (face.getDistanceAndNormal(*p, distance, normal,interactionRadius))
        {
            normal = -normal;
            BaseInteraction* const c = interactionHandler->getInteraction(p, this, timeStamp, normal);
            c->setNormal(normal);
            c->setDistance(distance);
            c->setOverlap(p->getRadius() - distance);
            c->setContactPoint(p->getPosition() - distance * normal);
            return c;
        }
    }
    return nullptr;
}

Mdouble TriangulatedWall::Face::getDistance(const Vec3D& otherPosition) const
{
    return Vec3D::dot(*vertex[0] - otherPosition, normal);
}

/**
 * check if there is contact with the face, determine if contact is with face, edge, vertex, return distance and normal;
 * only return edge, vertex contact if neighbor face pointer is higher to avoid doubles
 */
bool TriangulatedWall::Face::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return, Mdouble interactionRadius) const
{
    //check if particle overlaps
    distance = getDistance(p.getPosition());
    //return if distance is large, or particle has too much overlap
    if (fabs(distance) >= interactionRadius)///\todo make the triangle work from both sides
        return false;
    
    //Contact radius squared
    Mdouble allowedDistanceSquared = interactionRadius * interactionRadius - distance * distance;
    int touchingEdge = -1;
    Vec3D* touchingVertex = nullptr;
    //loop through edges to find out if there is a contact with face, edge, vertex, or neither.
    for (unsigned i = 0; i < 3; i++)
    {
        //distance from the edge
        Mdouble edgeDistance = Vec3D::dot(edgeNormal[i], *vertex[i] - p.getPosition());
        if (edgeDistance <= 0) //if in contact with the face
            continue;
        if (edgeDistance * edgeDistance >= allowedDistanceSquared) //if not in contact with anything
            return false;
        //this point is only reached if in contact, but not with the face
        //if edge is convex, treat as face, not edge contact
        bool convex = (neighbor[i] && ((neighbor[i]->getDistance(p.getPosition()) < 0) ? (
                Vec3D::dot(neighbor[i]->normal, edgeNormal[i]) > 1e-12) : (
                                               Vec3D::dot(neighbor[i]->normal, edgeNormal[i]) < -1e-12)));
        if (touchingEdge == -1)
        {
            //edge or vertex contact (depending if more edges are found)
            touchingEdge = i;
        }
        else
        {
            //vertex contact
            touchingVertex = vertex[(i == 2 && touchingEdge == 0) ? 0 : i];
        }
        //if (!convex) continue;
        //do not compute if neighbor exists and has lower index or face contact.
        if (neighbor[i] > this)
        { //if neighbor has higher index
            if (neighbor[i]->neighbor[0] == this)
            {
                edgeDistance = Vec3D::dot(neighbor[i]->edgeNormal[0], *vertex[i] - p.getPosition());
            }
            else if (neighbor[i]->neighbor[1] == this)
            {
                edgeDistance = Vec3D::dot(neighbor[i]->edgeNormal[1], *vertex[i] - p.getPosition());
            }
            else
            { //if (neighbor[i]->neighbor[2]==this)
                edgeDistance = Vec3D::dot(neighbor[i]->edgeNormal[2], *vertex[i] - p.getPosition());
            }
            if (edgeDistance <= 0 && !convex) //if neighbor has face contact, ignore the edge contact
                return false;
        }
        else if (neighbor[i] && !convex)
        {
            return false;
        }
    }
    
    
    //check if convex neighbours are overlapping as well
    if (touchingVertex)
    { //vertex contact
        normal_return = *touchingVertex - p.getPosition();
        distance = Vec3D::getLength(normal_return);
        normal_return /= distance;
        //return false;
        //logger(INFO,"vertex contact");
    }
    else if (touchingEdge == -1)
    { //face contact
        if (distance >= 0)
        { // front face contact
            normal_return = normal;
        }
        else
        { // back face contact
            normal_return = -normal;
            distance = -distance;
        }
        //return false;
        //logger(INFO,"face contact");
    }
    else
    { //edge contact
        Vec3D VP = *vertex[touchingEdge] - p.getPosition();
        Vec3D VW = *vertex[touchingEdge] - *vertex[(touchingEdge == 2) ? 0 : (touchingEdge + 1)];
        normal_return = VP - Vec3D::dot(VP, VW) / Vec3D::getLengthSquared(VW) * VW;
        distance = Vec3D::getLength(normal_return);
        normal_return /= distance;
        //return false;
        //logger(INFO,"edge contact at c=% p=% d=% n=%",p.getPosition() + distance * normal_return,p.getPosition(), distance, normal_return);
    }
    return true;
}

void TriangulatedWall::writeVTK(VTKContainer& vtk) const
{
    const int s = vtk.points.size();
    for (auto v : vertex_)
    {
        vtk.points.push_back(v);
    }
    for (auto f : face_)
    {
        std::vector<double> cell;
        cell.reserve(3);
        cell.push_back(s + f.vertex[0] - &vertex_[0]);
        cell.push_back(s + f.vertex[1] - &vertex_[0]);
        cell.push_back(s + f.vertex[2] - &vertex_[0]);
        vtk.triangleStrips.push_back(cell);
    }
}
