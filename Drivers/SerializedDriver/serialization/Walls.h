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

/* 
 * File:   Walls.h
 * Author: dducks
 *
 * Created on April 12, 2014, 2:19 PM
 */

#ifndef SERIALIZATION_WALLS_H
#define	SERIALIZATION_WALLS_H

#include <Walls/BaseWall.h>
#include <Walls/InfiniteWall.h>
#include <Walls/CylindricalWall.h>
#include <Walls/InfiniteWallWithHole.h>

namespace SerializationWrappers {
    template< typename Base >
    struct Wrapper {
        Wrapper() {
            data = nullptr;
        }
        ~Wrapper() {
            if (data != nullptr)
                delete data;
        }
        
        Base * data;
    };
}

template<class Archive>
void save(Archive& ar, const WallHandler& w) {
    ar ( cereal::make_size_tag(w.getNumberOfObjects()));
    for (const auto& wall : w ) {
        ar ( w );
    }
}

template<class Archive>
void load(Archive& ar, WallHandler& w) {
    cereal::size_type count;
    ar ( cereal::make_size_tag(count));
    
    //Due to inheritance etc etc, we need to create pointers here.
    std::cout << "WallHandler" << std::endl;
    for (int i = 0; i < count; i++) {
        SerializationWrappers::Wrapper<BaseWall> wall;
        std::cout << " Pass " << i << " / " << count << std::endl;
        ar ( wall );
        std::cout << " Adding..." << std::endl;
        w.copyAndAddObject( wall.data );
        std::cout << " Done. " << std::endl;
    }
}

template<class Archive>
void save(Archive& ar, const SerializationWrappers::Wrapper<BaseWall> w ) {
    if (typeid(*(w.data)) == typeid(InfiniteWall)) {
        ar( cereal::make_nvp("type","InfiniteWall"));
        save( ar, *(dynamic_cast<InfiniteWall*>(w.data)));
        //ar( cereal::make_nvp("value", dynamic_cast<const InfiniteWall*>(b)));
    /*} else if (typeid(b) == typeid(FiniteWall)) {
        ar( cereal::make_nvp("type","FiniteWall"));
        ar( cereal::make_nvp("value", dynamic_cast<const FiniteWall*>(b))); */
    } else if (typeid(*(w.data)) == typeid (InfiniteWallWithHole)) {
        ar( cereal::make_nvp("type","InfiniteWallWithHole"));
        //save( ar, *(dynamic_cast<InfiniteWallWithHole*>(w.data)));
    } else if (typeid(*(w.data)) == typeid (CylindricalWall)) {
        ar( cereal::make_nvp("type","CylindricalWall"));
        save( ar, *(dynamic_cast<CylindricalWall*>(w.data)));
    }
    
}

template<class Archive>
void load(Archive& ar, SerializationWrappers::Wrapper<BaseWall> & b ) {
    std::cout << "In load Generic Wall!" << std::endl;
    std::string type;
    ar( cereal::make_nvp("type", type));
    std::cout << "Type = " << type << std::endl;
    if (type == "InfiniteWall") {
        b.data = new InfiniteWall();
        //ar( cereal::make_nvp("value", dynamic_cast<InfiniteWall*>(b)));
        load( ar, *(dynamic_cast<InfiniteWall*>(b.data)));
   /* } else if (type == "FiniteWall") {
        b = new InfiniteWall();
        ar( cereal::make_nvp("value", dynamic_cast<FiniteWall*>(b)));*/
    } else if (type == "InfiniteWallWithHole") {
        b.data = new InfiniteWallWithHole();
        //ar( cereal::make_nvp("value", dynamic_cast<InfiniteWallWithHole*>(b)));
    } else if (type == "CylindricalWall") {
        b.data = new CylindricalWall();
        //ar( cereal::make_nvp("value", dynamic_cast<CylindricalWall*>(b)));
        load( ar, *(dynamic_cast<CylindricalWall*>(b.data)));
    }
    
}

template<class Archive>
void load(Archive& ar, InfiniteWall & w) {
    Vec3D position, normal;
    ar( CEREAL_NVP( position ),
        CEREAL_NVP( normal ));
    w.setPosition( position );
    w.setNormal( normal );
}

template<class Archive>
void save(Archive& ar, const InfiniteWall & w ) {
    ar( cereal::make_nvp("position", w.getPosition()),
        cereal::make_nvp("normal", w.getNormal()));
}

template<class Archive>
void load(Archive& ar, CylindricalWall & w ) {
    //Vec3D position;
    Mdouble radius;
    ar( //CEREAL_NVP( position ),
        CEREAL_NVP( radius ));
}

template<class Archive>
void save(Archive& ar, const CylindricalWall & w ) {
    ar ( //cereal::make_nvp("position", w.getPosition()),
         cereal::make_nvp("radius", w.getRadius()));
}

#endif	/* WALLS_H */

