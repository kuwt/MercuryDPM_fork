//
// Created by reza on 04-01-21.
//

#ifndef MERCURY_PIDCONTROLLER_H
#define MERCURY_PIDCONTROLLER_H
#include "Math/ExtendedMath.h"

class PIDController {
    Mdouble pGain;
    Mdouble iGain;
    Mdouble dGain;
    Mdouble iController = 0;
    bool filterFlag=true;
    Mdouble previousError;
    Mdouble filterError;

public:
// Constructor
    PIDController(Mdouble pGain, Mdouble iGain,Mdouble dGain)
            : pGain(pGain), iGain(iGain),dGain(iGain)
    {}
    Mdouble apply (Mdouble stressError, Mdouble timeStep);
};


#endif //MERCURY_PIDCONTROLLER_H
