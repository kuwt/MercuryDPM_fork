//
// Created by reza on 11-01-21.
//

#ifndef MERCURYDPM_PCONTROLLER_H
#define MERCURYDPM_PCONTROLLER_H
#include "Math/ExtendedMath.h"

class PController {
    Mdouble pGain;
public:

    PController(Mdouble pGain);

    Mdouble apply (Mdouble stressError, Mdouble timeStep);
};


#endif //MERCURYDPM_PCONTROLLER_H
