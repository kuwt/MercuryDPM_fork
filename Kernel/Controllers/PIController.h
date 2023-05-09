//
// Created by reza on 04-01-21.
//

#ifndef MERCURYDPM_PICONTROLLER_H
#define MERCURYDPM_PICONTROLLER_H
#include "Math/ExtendedMath.h"

class PIController {
    Mdouble pGain;
    Mdouble iGain;
    Mdouble iController = 0;

public:
// Constructor
    PIController(Mdouble pGain, Mdouble iGain);
    Mdouble  apply (Mdouble stressError, Mdouble timeStep);
};


#endif //MERCURYDPM_PICONTROLLER_H
