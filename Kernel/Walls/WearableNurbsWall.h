//Copyright (c) 2013-2021, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef MERCURYDPM_WEARABLENURBSWALL_H
#define MERCURYDPM_WEARABLENURBSWALL_H

#include "Walls/NurbsWall.h"

class WearableNurbsWall : public NurbsWall
{
public:
    
    /*!
     * \brief Default constructor: make a wall with default parameters.
     */
    WearableNurbsWall();
    
    /*!
     * \brief Copy constructor, copies another wall.
     */
    WearableNurbsWall(const WearableNurbsWall& other);
    
    /*!
     * \brief Constructor in which all parameters of the wall are set.
     */
    WearableNurbsWall(const NurbsSurface& nurbsSurface);
    
    /*!
     * \brief Default destructor.
     */
    ~WearableNurbsWall();
    
    WearableNurbsWall(Mdouble lengthU, Mdouble lengthV, Mdouble resolutionU, Mdouble resolutionV, bool periodicU = false, bool periodicV = false);
    
    void set(Mdouble lengthU, Mdouble lengthV, Mdouble resolutionU, Mdouble resolutionV, bool periodicU = false, bool periodicV = false);
    
    /*!
     * \brief Reads this wall from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes this wall to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object, here the string "WearableNurbsWall".
     */
    std::string getName() const final;
    
    /*!
     * \brief Copy this wall and return a pointer to the copy.
     */
    WearableNurbsWall* copy() const final;
    
    void computeWear() override;
    
    void writeWallDetailsVTK(VTKData& data) const override;
    
private:
    std::vector<std::vector<Mdouble>> localDebris_;
    void storeDebris(Vec3D P, const Mdouble debris);
    void processDebris();
public: // temp
    Mdouble getVolumeUnderSurface(const std::vector<Mdouble>& knotsU, const std::vector<Mdouble>& knotsV,
                                  const std::vector<std::vector<Vec3D>>& controlPoints, const std::vector<std::vector<Mdouble>>& weights) const;
    Mdouble getVolumeUnderSurfaceX(const std::vector<Mdouble>& knotsU, const std::vector<Mdouble>& knotsV,
                                   const std::vector<std::vector<Vec3D>>& controlPoints, const std::vector<std::vector<Mdouble>>& weights) const;
    //temp
public:
    void moveControlPoint(unsigned idxU, unsigned idxV, Vec3D dP)
    {
        nurbsSurface_.moveControlPoint(idxU, idxV, dP, false);
    }
    Mdouble getVolumeUnderSurface()
    {
        return getVolumeUnderSurface(nurbsSurface_.getKnotsU(), nurbsSurface_.getKnotsV(), nurbsSurface_.getControlPoints(), nurbsSurface_.getWeights());
    }
};


#endif //MERCURYDPM_WEARABLENURBSWALL_H
