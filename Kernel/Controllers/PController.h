//
// Created by reza on 11-01-21.
//

#ifndef MERCURY_PCONTROLLER_H
#define MERCURY_PCONTROLLER_H
#include "Math/ExtendedMath.h"

class PController {
    Mdouble pGain;
public:

    PController(Mdouble pGain);

    Mdouble apply (Mdouble stressError, Mdouble timeStep);
};


#endif //MERCURY_PCONTROLLER_H
