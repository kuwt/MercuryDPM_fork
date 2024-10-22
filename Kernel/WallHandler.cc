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


#include <Math/Helpers.h>
#include <Walls/BasicIntersectionOfWalls.h>
#include <Walls/TriangleWall.h>
#include <Walls/MeshTriangle.h>
#include <Walls/BasicUnionOfWalls.h>
#include <Walls/ScrewsymmetricIntersectionOfWalls.h>
#include "WallHandler.h"
#include "Walls/BaseWall.h"
#include "Walls/CylindricalWall.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/InfiniteWall.h"
#include "Walls/InfiniteWallWithHole.h"
#include "Walls/NurbsWall.h"
#include "Walls/WearableNurbsWall.h"
#include "Walls/Screw.h"
#include "Walls/Coil.h"
#include "DPMBase.h"
#include "Walls/VChute.h"
#include "BinaryReader.h"
#include "STLTriangle.h"
#include "Walls/TriangleMeshWall.h"
#include "Walls/WearableTriangleMeshWall.h"

/*!
 * Constructor of the WallHandler class. It creates an empty WallHandler.
 */
WallHandler::WallHandler()
{
    ///\todo why is this being done?
    clear();
    writeVTK_ = FileType::NO_FILE;
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "WallHandler::WallHandler() finished" << std::endl;
#endif
}

/*! 
 * \param[in] WH The WallHandler that has to be copied.
 * \details This is not a copy constructor! It only copies the pointer to the 
 *          DPMBase and the BaseWall in objects_, it sets the other data members
 *          to 0 or nullptr.
 */
WallHandler::WallHandler(const WallHandler& WH)
{
    ///\todo why is this being done?
    clear();
    setDPMBase(WH.getDPMBase());
    copyContentsFromOtherHandler(WH);
    writeVTK_ = WH.writeVTK_;
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "WallHandler::WallHandler(const WallHandler&) finished" << std::endl;
#endif
}

/*!
 * \param[in] rhs The WallHandler on the right hand side of the assignment.
 * \details This is not a copy assignment operator! It only copies the pointer to the 
 *          DPMBase and the BaseWall in objects_, it sets the other data members
 *          to 0 or nullptr.
 */
WallHandler& WallHandler::operator=(const WallHandler& rhs)
{
    if (this != &rhs)
    {
        clear();
        setDPMBase(rhs.getDPMBase());
        copyContentsFromOtherHandler(rhs);
        writeVTK_ = rhs.writeVTK_;
    }
    return *this;
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "WallHandler::operator =(const WallHandler&) finished" << std::endl;
#endif
}

WallHandler::~WallHandler()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "WallHandler::~WallHandler() finished" << std::endl;
#endif
}

/*!
 * \param[in] W A pointer to the BaseWall (or derived class) that has to be added.
 * \details First the new BaseWall is added to the vector of BaseWall, then it is 
 * told that this is its handler.
 */
void WallHandler::addObject(BaseWall* W)
{
    if (W->getSpecies() == nullptr)
    {
        logger(WARN, "WARNING: The wall with ID % that is added in WallHandler::addObject "
                     "does not have a species yet. Please make sure that you have "
                     "set the species somewhere in the driver code.", W->getId());
    }
    //Puts the wall in the Wall list
    BaseHandler<BaseWall>::addObject(W);
    //set the particleHandler pointer
    W->setHandler(this);
}

BaseWall* WallHandler::createObject(const std::string& type)
{
    if (type == "CylindricalWall")
    {
        return new CylindricalWall;
    }
    else if (type == "AxisymmetricIntersectionOfWalls")
    {
        return new AxisymmetricIntersectionOfWalls;
    }
    else if (type == "ScrewsymmetricIntersectionOfWalls")
    {
        return new ScrewsymmetricIntersectionOfWalls;
    }
    else if (type == "IntersectionOfWalls")
    {
        return new IntersectionOfWalls;
    }
    else if (type == "BasicIntersectionOfWalls")
    {
        return new BasicIntersectionOfWalls;
    }
    else if (type == "BasicUnionOfWalls")
    {
        return new BasicUnionOfWalls;
    }
    else if (type == "InfiniteWall")
    {
        return new InfiniteWall;
    }
    else if (type == "InfiniteWallWithHole")
    {
        return new InfiniteWallWithHole;
    }
    else if (type == "Screw")
    {
        return new Screw;
    }
    else if (type == "Coil")
    {
        return new Coil;
    }
    else if (type == "TriangleWall")
    {
        return new TriangleWall;
    }
    else if (type == "MeshTriangle")
    {
        return new MeshTriangle;
    }
    else if (type == "VChute")
    {
        return new VChute;
    }
    else if (type == "NurbsWall")
    {
        return new NurbsWall();
    }
    else if (type == "WearableNurbsWall")
    {
        return new WearableNurbsWall();
    }
    else if (type == "TriangleMeshWall")
    {
        return new TriangleMeshWall;
    }
    else if (type == "WearableTriangleMeshWall")
    {
        return new WearableTriangleMeshWall;
    }
    //for backward compatibility (before svnversion ~2360)
    else if (type == "numFiniteWalls")
    {
        return new BasicIntersectionOfWalls;
    }
        /// \todo Review this line. Problem came up in merging.
   // else if (!getDPMBase()->readUserDefinedWall(type,is))
    else
    {
        logger(WARN, "Wall type: % not understood in restart file", type);
        return nullptr;
    }
}

BaseWall* WallHandler::readAndCreateObject(std::istream& is)
{
    std::string type;
    is >> type;
    logger(DEBUG, "WallHandler::readAndAddObject(is): reading type %.", type);

    //for backward compatibility (before svnversion ~2360)
    if (type == "numFiniteWalls")
    {
        return readAndCreateOldObject(is);
    }
    else
    {
        BaseWall* wall = createObject(type);
        //check if wall is user-defined
        if (wall == nullptr)
        {
            wall = getDPMBase()->readUserDefinedWall(type);
        }
        //throw warning if wall could not be found
        if (wall == nullptr)
        {
            std::string line;
            getline(is, line);
            logger(WARN, "This wall could not be read; dummy wall is inserted instead:\n%%", type, line);
            BaseWall* wall = new InfiniteWall;
            wall->setHandler(this);
            wall->setSpecies(getDPMBase()->speciesHandler.getObject(0));
            return wall;
        }
        wall->setHandler(this);
        is >> *wall;
        wall->setSpecies(getDPMBase()->speciesHandler.getObject(wall->getIndSpecies()));
        return wall;
    }
}

/*!
 * \details First determine whether or not the wall is an infinite wall. If it is
 * an infinite wall, read the normal and position and add the wall to the handler.
 * If it is a  finite wall, read the normal and position of each part and construct
 * an IntersectionOfWalls from it, which can then be added to the handler.
 * \param[in,out] is The input stream from which the information is read.
 *
 * \todo This is deprecated since r ~2360.
 */
BaseWall* WallHandler::readAndCreateOldObject(std::istream& is)
{
    //read in next line
    std::stringstream line;
    helpers::getLineFromStringStream(is, line);
    logger(VERBOSE, line.str());

    std::string dummy;
    unsigned int numWalls;
    Mdouble position;
    Vec3D normal;
    line >> numWalls;

    if (numWalls == 0)
    {
        InfiniteWall* wall = new InfiniteWall();
        wall->setSpecies(getDPMBase()->speciesHandler.getObject(0));
        line >> dummy >> normal >> dummy >> position;
        wall->set(normal, position * normal);
        return wall;
    }
    else
    {
        IntersectionOfWalls* wall = new IntersectionOfWalls();
        wall->setSpecies(getDPMBase()->speciesHandler.getObject(0));
        for (unsigned int i = 0; i < numWalls; ++i)
        {
            line >> dummy >> normal >> dummy >> position;
            wall->addObject(normal, position * normal);
        }
        return wall;
    }
}

/*!
 * \details As we add the object into the handler, we need to make sure that the
 * object keeps its existing ID. This is not always the same as the index, e.g.
 * if walls have been removed during the simulation.
 * \param[in] is The input stream from which the information is read.
 */
void WallHandler::readAndAddObject(std::istream& is)
{
    BaseWall* o = readAndCreateObject(is);
    unsigned int id = o->getId();
    addObject(o);
    getLastObject()->setId(id);
}

/*!
 * \return The string "WallHandler".
 */
std::string WallHandler::getName() const
{
    return "WallHandler";
}

void WallHandler::writeWallDetailsVTKBoundingBox(VTKData& data) const
{
    Vec3D P[2] = { getDPMBase()->getMax(), getDPMBase()->getMin() };
    for (const Vec3D& i : P)
    {
        for (const Vec3D& j : P)
        {
            for (const Vec3D& k : P)
            {
                data.addToPoints(Vec3D(i.X, j.Y, k.Z));
            }
        }
    }

    data.addToConnectivity({ 0, 1, 3, 2, 0, 4, 5, 1, 5, 7, 3, 7, 6, 2, 6, 4 });
    data.addToTypes(4);
}

/**
 * \details The stl files have to be binary STL files. If you have ascii files, you need to convert them (use e.g. https://www.meshconvert.com/)
 *
 * The vtk files need to be of type POLYDATA and contain triangle strips (see <a href="http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html">see www.cacr.caltech.edu/~slombey</a>)
 * \param[in] filename name of vtk input file, e.g. TriangulatedWallSelfTest.vtk
 * \param[in] species pointer to a species in the species handler that will be assigned to the walls
 * \param[in] scaleFactor allows the vertex positions to be scaled (necessary if the vtk file is written in different units than the Mercury implementation, e.g. if the stl file is given in mm, but the Mercury implementation uses meters)
 */
unsigned WallHandler::readTriangleWall(std::string filename, ParticleSpecies* species, Mdouble scaleFactor, Vec3D centerOfRotation, Vec3D velocity, Vec3D angularVelocity)
{
    const unsigned groupId = getNextGroupId();
    std::string fileType = filename.substr(filename.find_last_of('.') + 1);

    //define a default triangle wall
    TriangleWall triangleWall;
    triangleWall.setSpecies(species);
    triangleWall.setVelocity(velocity);
    triangleWall.setAngularVelocity(angularVelocity);
    triangleWall.setGroupId(groupId);

    if (helpers::lower(fileType) == "vtk")
    {
        //try open the input file
        std::fstream file;
        file.open(filename.c_str(), std::ios::in);
        logger.assert_always(file.is_open(), "File opening failed: %", filename);

        //skip the header lines
        std::string dummy;
        getline(file, dummy);
        getline(file, dummy);
        getline(file, dummy);
        getline(file, dummy);

        //read vertices, apply scaling
        unsigned num;
        file >> dummy >> num >> dummy;
        std::vector<Vec3D> vertex;
        vertex.reserve(num);
        Vec3D v;
        for (unsigned i = 0; i < num; i++)
        {
            file >> v.X >> v.Y >> v.Z;
            v *= scaleFactor;
            vertex.push_back(v);
        }

        //read faces
        unsigned n = getSize();
        file >> dummy >> num >> dummy;
        unsigned id0, id1, id2;
        for (unsigned i = 0; i < num; i++)
        {
            file >> dummy >> id0 >> id1 >> id2;
            triangleWall.setVertices(vertex[id0], vertex[id1], vertex[id2], centerOfRotation);
            copyAndAddObject(triangleWall);
        }

        //close file
        file.close();

        logger(INFO, "Read in % walls from %", getSize() - n,filename);

    }
    else if (helpers::lower(fileType) == "stl")
    {

        BinaryReader file(filename);

        STLTriangle triangle;

        std::string header = file.readString(80);
        unsigned numTriangles = file.readUnsignedInt(4);

        for (unsigned i = 0; i < numTriangles; i++)
        {
            triangle.normal.x() = file.readFloat(4);
            triangle.normal.y() = file.readFloat(4);
            triangle.normal.z() = file.readFloat(4);


            triangle.vertex1.x() = file.readFloat(4);
            triangle.vertex1.y() = file.readFloat(4);
            triangle.vertex1.z() = file.readFloat(4);

            triangle.vertex2.x() = file.readFloat(4);
            triangle.vertex2.y() = file.readFloat(4);
            triangle.vertex2.z() = file.readFloat(4);


            triangle.vertex3.x() = file.readFloat(4);
            triangle.vertex3.y() = file.readFloat(4);
            triangle.vertex3.z() = file.readFloat(4);

            triangle.vertex1 *= scaleFactor;
            triangle.vertex2 *= scaleFactor;
            triangle.vertex3 *= scaleFactor;

            //add to triangle wall
            triangleWall.setVertices(triangle.vertex1, triangle.vertex2, triangle.vertex3, centerOfRotation);
            copyAndAddObject(triangleWall);

            //Now ignore (read) the two dummy characters
            file.ignoreChar(2);

        }

        logger(INFO, "Read in % walls from %", numTriangles,filename);

    }
    else
    {

        logger(ERROR, "File type of % must be vtk or stl");

    }

    return groupId;
}

/**
 * \details This function is called by DPMBase::computeOneTimeStep directly after
 * the ghost particle update. It may be used to adjust wall properties based on
 * the changed wall positions (for an example, look e.g. at MeshTriangle).
 */
void WallHandler::actionsAfterParticleGhostUpdate()
{
    // #pragma omp parallel for //schedule(dynamic)
    for (auto w: *this)
    {
        w->actionsAfterParticleGhostUpdate();
    }
}

void WallHandler::setWriteVTK(FileType fileType)
{
    writeVTK_ = fileType;
}

void WallHandler::setWriteVTK(bool writeVTK)
{
    writeVTK_ = writeVTK ? FileType::MULTIPLE_FILES : FileType::NO_FILE;
}

FileType WallHandler::getWriteVTK() const
{
    return writeVTK_;
}

void WallHandler::setWriteDetailsVTK(DetailsVTKOptions option, FileType fileType)
{
    // Automatically adds option when it doesn't exist in the map yet
    writeDetailsVTK_[option] = fileType;
}

void WallHandler::setWriteDetailsVTK(DetailsVTKOptions option, bool writeVTK)
{
    FileType fileType = writeVTK ? FileType::MULTIPLE_FILES : FileType::NO_FILE;
    setWriteDetailsVTK(option, fileType);
}

FileType WallHandler::getWriteDetailsVTK(DetailsVTKOptions option) const
{
    // This function has to be marked constant, so accessing by key is not an option, as it will automatically add keys
    // not present in the map yet.

    // Try to find the option to see if it's already added to the map or not
    auto it = writeDetailsVTK_.find(option);
    // When the iterator does not equal the end, the key has been found
    if (it != writeDetailsVTK_.end())
        return it->second;
    else
        return FileType::NO_FILE; // Else return a default file type
}

bool WallHandler::getWriteDetailsVTKAny() const
{
    // Simple lambda expression to check if any of the file types are not of type NO_FILE
    return std::any_of(writeDetailsVTK_.begin(), writeDetailsVTK_.end(),
                       [](const auto& p) { return p.second != FileType::NO_FILE; });
}

std::unordered_map<WallHandler::DetailsVTKOptions, FileType> WallHandler::getWriteWallDetailsVTKAll() const
{
    return writeDetailsVTK_;
}

void WallHandler::setWriteWallSurfaceAreaVTK(bool writeWallSurfaceAreaVTK)
{
    getDPMBase()->getWallVTKWriter().setWriteWallSurfaceAreaVTK(writeWallSurfaceAreaVTK);
}

bool WallHandler::getWriteWallSurfaceAreaVTK() const 
{
    return getDPMBase()->getWallVTKWriter().getWriteWallSurfaceAreaVTK();
}
