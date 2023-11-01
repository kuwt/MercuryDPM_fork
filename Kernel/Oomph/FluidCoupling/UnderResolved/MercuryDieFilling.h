//
// Created by mitchel on 6/4/19.
//  Modified by Hao on 2-Jul-2019
//

#ifndef MERCURYDPM_MERCURYDIEFILLING_H
#define MERCURYDPM_MERCURYDIEFILLING_H

// Generic MercuryDPM header
#include "Mercury3D.h"
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include <Walls/InfiniteWall.h>

#include "../../../../oomph-lib/src/generic/generic.h"

#include <cmath>

//namespace mercVoidage
//{
    //Mercury3D *dpmPointer;
    //double getVoidageFromMerc(const double &, const oomph::Vector<double> &);
//}

class MercuryDieFilling : public Mercury3D
{
public:
    MercuryDieFilling() {}
    ~MercuryDieFilling() {}
};

#endif //MERCURYDPM_MERCURYDIEFILLING_H