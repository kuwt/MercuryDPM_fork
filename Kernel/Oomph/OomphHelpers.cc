#include "OomphHelpers.h"

namespace convertVecFuncs
{
    Vec3D convertToVec3D(oomph::Vector<double>& Vec)
    {
        Vec3D returnVec;
        returnVec.X = Vec[0];
        returnVec.Y = Vec[1];
        returnVec.Z = Vec[2];
        return returnVec;
    }

    oomph::Vector<double> convertToOomphVec(const Vec3D& Vec)
    {
        oomph::Vector<double> returnVec(3,0.0);
        returnVec[0] = Vec.X;
        returnVec[1] = Vec.Y;
        returnVec[2] = Vec.Z;
        return returnVec;
    }
}