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

#ifndef MERCURY_WEARABLETRIANGULATEDWALL_H
#define MERCURY_WEARABLETRIANGULATEDWALL_H

#include "TriangulatedWall.h"

class WearableTriangulatedWall : public TriangulatedWall
{
public:
    
    /*!
     * \brief Default constructor.
     */
    WearableTriangulatedWall();
    
    /*!
     * \brief Copy constructor.
     */
    WearableTriangulatedWall(const WearableTriangulatedWall& other);
    
    /*!
     * \brief Creates flat triangulated rectangle with certain resolution.
     * @param lengthU Length in u-direction
     * @param lengthV Length in v-direction
     * @param resolutionU Max stepsize in u-direction
     * @param resolutionV Max stepsize in v-direction
     */
    WearableTriangulatedWall(Mdouble lengthU, Mdouble lengthV, Mdouble resolutionU, Mdouble resolutionV);
    
    /*!
     * \brief Destructor.
     */
    ~WearableTriangulatedWall() override;
    
    /*!
     * Copy assignment operator.
     */
    WearableTriangulatedWall& operator=(const WearableTriangulatedWall& other);
    
    /*!
     * \brief Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
     */
    WearableTriangulatedWall* copy() const override;
    
    /*!
     * \brief Reads an WearableTriangulatedWall from an input stream, for example a restart file.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes an WearableTriangulatedWall to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object, here the string "WearableTriangulatedWall".
     */
    std::string getName() const override;
    
    void computeWear() override;

private:
    std::vector<std::vector<Vec3D>> debris_;

    void storeDebris(Vec3D P, const Vec3D& debris);
    
    void processDebris();
};


#endif //MERCURY_WEARABLETRIANGULATEDWALL_H
