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
 * File:   Problem.h
 * Author: dducks
 *
 * Created on April 12, 2014, 2:20 PM
 */

#ifndef SERIALIZATION_PROBLEM_H
#define	SERIALIZATION_PROBLEM_H

#include <cereal/cereal.hpp>
#include <Mercury3D.h>

class SerializedProblem : public Mercury3D {
public:
    //SerializedProblem();
};

template<class Archive>
void load(Archive& ar, SerializedProblem & base) {
    using namespace cereal;
    
    std::string name;
    Mdouble tMax;
    Mdouble tStep;
    Vec3D gravity;
    Vec3D minimum;
    Vec3D maximum;
        
    ar( CEREAL_NVP( name )
        ,CEREAL_NVP( tMax )
        ,CEREAL_NVP( tStep )
        ,CEREAL_NVP( gravity )
        ,CEREAL_NVP( minimum )
        ,CEREAL_NVP( maximum )
        ,cereal::make_nvp("particles", base.particleHandler)
        ,cereal::make_nvp("walls", base.wallHandler)
    );
    
    base.setName( name );
    base.setTimeMax( tMax );
    base.setTimeStep( tStep );
    base.setGravity( gravity );
    base.setXMin(minimum.X);
    base.setYMin(minimum.Y);
    base.setZMin(minimum.Z);
    
    base.setXMax(maximum.X);
    base.setYMax(maximum.Y);
    base.setZMax(maximum.Z);
}

template<class Archive>
void save(Archive& ar, const SerializedProblem & base) {
    Vec3D minimum;
    Vec3D maximum;
    
    minimum.X = base.getXMin();
    minimum.Y = base.getYMin();
    minimum.Z = base.getZMin();
    maximum.X = base.getXMax();
    maximum.Y = base.getYMax();
    maximum.Z = base.getZMax();
    
    ar ( cereal::make_nvp("name", base.getName())
         ,CEREAL_NVP(minimum)
         ,CEREAL_NVP(maximum)
         ,cereal::make_nvp("tMax", base.getTimeMax())
         ,cereal::make_nvp("tStep", base.getTimeStep())
         ,cereal::make_nvp("gravity", base.getGravity())
         ,cereal::make_nvp("particles", base.particleHandler)
         ,cereal::make_nvp("walls", base.wallHandler)
       );
}


#endif	/* PROBLEM_H */

