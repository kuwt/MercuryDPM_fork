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


#ifndef WALLHANDLER_H
#define WALLHANDLER_H

#include "BaseHandler.h"
#include "Walls/BaseWall.h"

class BaseWall;

/*!
 * \class WallHandler
 * \brief Container to store all BaseWall
 * 
 * \details The WallHandler is a container to store all BaseWall. 
 * It is implemented by a vector of pointers to BaseWall.
 */
class WallHandler final : public BaseHandler<BaseWall>
{
public:
    /*!
     * \brief Default constructor, it creates an empty WallHandler.
     */
    WallHandler();

    /*!
     * \brief Constructor that copies a pointer to the DPMBase and all BaseWall in the given WallHandler.
     */
    WallHandler(const WallHandler& BH);

    /*!
     * \brief Assignment operator that copies the pointer to the DPMBase and all BaseWall in the given WallHandler.
     */
    WallHandler& operator=(const WallHandler& rhs);
    
    /*!
     * \brief Destructor, it destructs the WallHandler and all BaseWall it contains.
     */
    ~WallHandler() final;

    /*!
     * \brief Adds a BaseWall to the WallHandler.
     */
    void addObject(BaseWall* W) final;

    /*!
     * \brief Create a new wall, with the type given as a string (required for restarting).
     */
    static BaseWall* createObject(const std::string& type);

    /*!
     * \brief Create a new wall, based on the information provided in a restart file.
     */
    BaseWall* readAndCreateObject(std::istream& is);

    /*!
     * \brief Create a new wall in the WallHandler, based on the information provided in a restart file.
     */
    void readAndAddObject(std::istream& is) final;

    /*!
     * \brief Returns the name of the handler, namely the string "WallHandler".
     */
    std::string getName() const final;

    /**
     * \brief Reads triangulated walls from vtk or stl files, and converts them into a set of TriangleWalls
     */

    unsigned readTriangleWall(std::string filename, ParticleSpecies* species, Mdouble scaleFactor = 1, Vec3D centerOfRotation = {0,0,0}, Vec3D velocity = {0,0,0}, Vec3D angularVelocity = {0,0,0});

private:

    /*!
     * \brief Create a new wall, based on the information from old-style restart data.
     */
    BaseWall* readAndCreateOldObject(std::istream& is);

    /*!
     * \brief Writes a bounding box around the domain into a vtk file.
     */
    void writeVTKBoundingBox() const;

};

#endif

