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

#include "Mercury3D.h"
#include "StatisticsVector.h"
#include "Walls/InfiniteWall.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <Species/LinearViscoelasticSpecies.h>

/// In this file a cubic packing of 5^3 particles in a tri-axial box is created and allowed to settle under small gravity. After that Z statistics are calculated.
class DPM : public Mercury3D
{

public:

    void setupInitialConditions() override {
        const double d = 1.;
        const unsigned N = 5;

        setMin({0,0,0});
        setMax({N*d,N*d,N*d});
        setTimeStep(1);
        speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        setName("ComputeVolumeFractionSelfTest");

        SphericalParticle P0;
        P0.setSpecies(speciesHandler.getObject(0));
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                for (int k = 0; k < N; k++)
                {
                    P0.setRadius(0.5*d);
                    P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
                    P0.setPosition(d * Vec3D(.5 + i, .5 + j, .5 + k));
                    particleHandler.copyAndAddObject(P0);
                }
    }

    /**
     * Computes the local volume fraction at a given position, using a Heaviside CG function with a given cutoff
     */
    Mdouble computeLocalVolumeFraction(Vec3D position, Mdouble cutoff) {
        const Mdouble widthSquared = cutoff*cutoff;
        const Mdouble inverseVolume = 1./(4./3.*constants::pi*widthSquared*cutoff);
        Mdouble volumeFraction = 0;
        for (auto p: particleHandler) {
            const Mdouble radialDistanceSquared = Vec3D::getLengthSquared(p->getPosition()-position);
            const Mdouble cgFunction = (radialDistanceSquared<widthSquared)*inverseVolume;
            volumeFraction += p->getVolume()*cgFunction;
        }
        return volumeFraction;
    }

    Mdouble computeLocalVolumeFractionHGrid(Vec3D position, Mdouble cutoff) {
        const Mdouble widthSquared = cutoff*cutoff;
        const Mdouble inverseVolume = 1./(4./3.*constants::pi*widthSquared*cutoff);
        Mdouble volumeFraction = 0;

        const Vec3D min = position - Vec3D(cutoff,cutoff,cutoff);
        const Vec3D max = position + Vec3D(cutoff,cutoff,cutoff);
        HGrid* const hGrid = getHGrid();
        unsigned int level = 0;
        for (int occupiedLevelsMask = hGrid->getOccupiedLevelsMask(); occupiedLevelsMask!=0; occupiedLevelsMask >>= 1, ++level)
        {
            // If no objects at this level, go on to the next level
            if ((occupiedLevelsMask & 1) == 0) continue;

            const Mdouble inv_size = hGrid->getInvCellSize(level);
            const int xs = static_cast<int>(std::floor(min.X * inv_size - 0.5));
            const int xe = static_cast<int>(std::floor(max.X * inv_size + 0.5));
            const int ys = static_cast<int>(std::floor(min.Y * inv_size - 0.5));
            const int ye = static_cast<int>(std::floor(max.Y * inv_size + 0.5));
            const int zs = static_cast<int>(std::floor(min.Z * inv_size - 0.5));
            const int ze = static_cast<int>(std::floor(max.Z * inv_size + 0.5));

            for (int x = xs; x <= xe; ++x)
            {
                for (int y = ys; y <= ye; ++y)
                {
                    for (int z = zs; z <= ze; ++z)
                    {
                        // Loop through all objects in the bucket to find nearby objects
                        const unsigned int bucket = hGrid->computeHashBucketIndex(x, y, z, level);
                        BaseParticle* p = hGrid->getFirstBaseParticleInBucket(bucket);
                        while (p != nullptr && p->getHGridCell().equals(x, y, z, level))
                        {
                            const Mdouble radialDistanceSquared = Vec3D::getLengthSquared(p->getPosition()-position);
                            const Mdouble cgFunction = (radialDistanceSquared<widthSquared)*inverseVolume;
                            volumeFraction += p->getVolume()*cgFunction;
                            p = p->getHGridNextObject();
                        }
                    }
                }
            }
        } //end for level

        return volumeFraction;
    }

    template <class T>
    Mdouble computeLocalCGHGrid(Vec3D position, Mdouble cutoff, T cgValue(BaseParticle*)) {
        const Mdouble widthSquared = cutoff*cutoff;
        const Mdouble inverseVolume = 1./(4./3.*constants::pi*widthSquared*cutoff);
        T value = 0;

        const Vec3D min = position - Vec3D(cutoff,cutoff,cutoff);
        const Vec3D max = position + Vec3D(cutoff,cutoff,cutoff);
        HGrid* const hGrid = getHGrid();
        unsigned int level = 0;
        for (int occupiedLevelsMask = hGrid->getOccupiedLevelsMask(); occupiedLevelsMask!=0; occupiedLevelsMask >>= 1, ++level)
        {
            // If no objects at this level, go on to the next level
            if ((occupiedLevelsMask & 1) == 0) continue;

            const Mdouble inv_size = hGrid->getInvCellSize(level);
            const int xs = static_cast<int>(std::floor(min.X * inv_size - 0.5));
            const int xe = static_cast<int>(std::floor(max.X * inv_size + 0.5));
            const int ys = static_cast<int>(std::floor(min.Y * inv_size - 0.5));
            const int ye = static_cast<int>(std::floor(max.Y * inv_size + 0.5));
            const int zs = static_cast<int>(std::floor(min.Z * inv_size - 0.5));
            const int ze = static_cast<int>(std::floor(max.Z * inv_size + 0.5));

            for (int x = xs; x <= xe; ++x)
            {
                for (int y = ys; y <= ye; ++y)
                {
                    for (int z = zs; z <= ze; ++z)
                    {
                        // Loop through all objects in the bucket to find nearby objects
                        const unsigned int bucket = hGrid->computeHashBucketIndex(x, y, z, level);
                        BaseParticle* p = hGrid->getFirstBaseParticleInBucket(bucket);
                        while (p != nullptr && p->getHGridCell().equals(x, y, z, level))
                        {
                            const Mdouble radialDistanceSquared = Vec3D::getLengthSquared(p->getPosition()-position);
                            const Mdouble cgFunction = (radialDistanceSquared<widthSquared)*inverseVolume;
                            value += cgValue(p)*cgFunction;
                            p = p->getHGridNextObject();
                        }
                    }
                }
            }
        } //end for level

        return value;
    }

    Mdouble computeLocalVolumeFractionHGridCompact(Vec3D position, Mdouble cutoff) {
        auto getVolume = [] (BaseParticle* p) {return p->getVolume();};
        return computeLocalCGHGrid<double> (position, cutoff, getVolume);
    }

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    DPM problem;
    problem.solve();
    Vec3D position = 2.5*Vec3D(1,1,1);
    double volumeFraction = problem.computeLocalVolumeFraction(position,2.5);
    logger(INFO,"LocalVolumeFraction at %: %",position, volumeFraction);
    volumeFraction = problem.computeLocalVolumeFractionHGrid(position,2.5);
    logger(INFO,"LocalVolumeFraction with HGrid at %: %",position, volumeFraction);
    volumeFraction = problem.computeLocalVolumeFractionHGridCompact(position,2.5);
    logger(INFO,"LocalVolumeFraction with HGrid at %: %",position, volumeFraction);
}

