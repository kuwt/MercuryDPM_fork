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

namespace mercVoidage
{
    Mercury3D *dpmPointer = nullptr;
    
    // USED FOR SEMI-RESOLVED ONLY
    double getVoidageFromMerc(const double &time, const oomph::Vector<double> &x)
    {
        if (time != 0 && time!= 1e-8)
        {
            std::cout << "time takes value it should not take: " << time << std::endl;
        }
        
        unsigned int nP = dpmPointer->particleHandler.getNumberOfObjects();
        
        
        double localVoidage = 1.0;
        //std::cout << "nP = " << nP << std::endl;
        
        for (unsigned iP = 0; iP < nP; iP++)
        {
            // Voidage at a certain position in time is equal to the voidage at that position - velocity*dt
            // Hence for get_voidage_nst and get_voidage_gradient_nst the dt term is zero, as function calls use equal time
            // For get_voidage_gradient_time_nst however time is different as t_pls = t+epsfd
            // time > dpmPointer->getTime(), thus if velocity is positive, currX is on the left of getPosition().X
            
            //FIXME time does not equal the global time of the simulation, but is equal to epsfd
            // hence *(time- dpmPointer->getTime()) is replaced by *(time).
            // time equals either 0 or epsfd, time is only for get_voidage_gradient_time_nst non-zero
            // note that time() is used as global time, but that can not be accessed inside this namespace
            
            double currX = dpmPointer->particleHandler.getObject(iP)->getPosition().X +
                           dpmPointer->particleHandler.getObject(iP)->getVelocity().X * (time);
            double currY = dpmPointer->particleHandler.getObject(iP)->getPosition().Y +
                           dpmPointer->particleHandler.getObject(iP)->getVelocity().Y * (time);
            double currZ = dpmPointer->particleHandler.getObject(iP)->getPosition().Z +
                           dpmPointer->particleHandler.getObject(iP)->getVelocity().Z * (time);
            
            double dist = pow((pow(x[0] - currX, 2) + pow(x[1] - currY, 2)), 0.5); // FIXME for 3D
            double c = 3.0;
            double cutoff = dpmPointer->particleHandler.getObject(iP)->getRadius();
            
            // Do note the fudge factor, as the code crashes for large particles
            localVoidage += 0.1 +  0.9 * ( 0.5 * (erf(c * (dist - cutoff))) - 0.5 * (erf(c * (dist + cutoff))));
        }
        return localVoidage;
    }
}


class MercuryDieFilling : public Mercury3D //, GeomObject
// Can not place GeomObject simply as it needs a function definition for
// oomph::GeomObject::position(const oomph::Vector<double>&, oomph::Vector<double>)
// which is declared pure virtual
{
public:
    MercuryDieFilling()
    {
        mercVoidage::dpmPointer = this;
    }
    ~MercuryDieFilling() {}
protected:

};

#endif //MERCURYDPM_MERCURYDIEFILLING_H
