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

#ifndef MERCURY_HGRIDCELL_H
#define MERCURY_HGRIDCELL_H

/// Contains the hGrid-information for a certain particle: x,y,z and level of the particle containing this.
/// Note, that each particle contains a HGridCell.
/// All methods are inline for speed reasons: they are called VERY often, so we do need performance here.
class HGridCell
{
public:
    HGridCell() : hGridX_(0), hGridY_(0), hGridZ_(0), hGridLevel_(0)
    {}
    
    ///Checks if the given (x,y,z,level) is the same as the ones in this cell.
    inline bool equals(int x, int y, int z, unsigned int level) const
    {
        return (x == hGridX_ && y == hGridY_ && z == hGridZ_ && level == hGridLevel_);
    }
    
    ///Checks if the given (x,y,z,level) is the same as the ones in this cell, 2D version.
    inline bool equals(int x, int y, unsigned int level) const
    {
        return (x == hGridX_ && y == hGridY_ && level == hGridLevel_);
    }
    
    ///Checks if the given cell is the same as the given cell
    inline bool operator==(const HGridCell& other) const
    {
        return equals(other.hGridX_, other.hGridY_, other.hGridZ_, other.hGridLevel_);
    }
    
    inline int getHGridX() const
    {
        return hGridX_;
    }
    
    inline void setHGridX(int HGridX)
    {
        hGridX_ = HGridX;
    }
    
    inline int getHGridY() const
    {
        return hGridY_;
    }
    
    inline void setHGridY(int HGridY)
    {
        hGridY_ = HGridY;
    }
    
    inline int getHGridZ() const
    {
        return hGridZ_;
    }
    
    inline void setHGridZ(int HGridZ)
    {
        hGridZ_ = HGridZ;
    }
    
    inline unsigned int getHGridLevel() const
    {
        return hGridLevel_;
    }
    
    inline void setHGridLevel(unsigned int HGridLevel)
    {
        hGridLevel_ = HGridLevel;
    }

private:
    
    ///Cell position in the grid
    int hGridX_, hGridY_, hGridZ_;
    ///HGrid-level of the particle containing this cell
    unsigned int hGridLevel_;
    
};


#endif //MERCURY_HGRIDCELL_H
