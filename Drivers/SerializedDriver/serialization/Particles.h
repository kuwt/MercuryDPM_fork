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
 * File:   Particles.h
 * Author: dducks
 *
 * Created on April 12, 2014, 2:19 PM
 */

#ifndef SERIALIZATION_PARTICLES_H
#define	SERIALIZATION_PARTICLES_H

#include <cereal/cereal.hpp>

#include <ParticleHandler.h>
#include <Particles/BaseParticle.h>

template<class Archive>
void load(Archive& ar, ParticleHandler& handl) {
    cereal::size_type size;
    ar ( cereal::make_size_tag(size) );
    
    SphericalParticle p;
    
    for (int i = 0; i < size; i++) {
        ar ( p );
        std::cout << "Particle:"
                  << "\n\tPOS: " << p.getPosition()
                  << "\n\tVEL: " << p.getVelocity()
                  << "\n\tRAD: " << p.getRadius() << std::endl;
        handl.copyAndAddObject( p );
    }
}

template<class Archive>
void save(Archive& ar, const ParticleHandler& handl) {
    ar ( cereal::make_size_tag( handl.getNumberOfObjects() ));
    for (const auto& p : handl) {
        ar ( p );
    }
}

template<class Archive>
void load(Archive& ar, BaseParticle& p) {
    Vec3D position;
    Vec3D velocity;
    Mdouble radius;
    
    ar( CEREAL_NVP(position),
        CEREAL_NVP(velocity),
        CEREAL_NVP(radius));
    
    p.setPosition(position);
    p.setRadius(radius);
    p.setVelocity(velocity);
}

template<class Archive>
void save(Archive& ar, const BaseParticle& p) {
    ar( cereal::make_nvp("position", p.getPosition()),
        cereal::make_nvp("velocity", p.getVelocity()),
        cereal::make_nvp("radius", p.getRadius()));
}



#endif	/* PARTICLES_H */

