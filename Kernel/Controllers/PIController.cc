//
// Created by reza on 04-01-21.
//

#include "PIController.h"
PIController::PIController(Mdouble pGain, Mdouble iGain)
: pGain(pGain), iGain(iGain)
{}
Mdouble PIController::apply (Mdouble stressError, Mdouble timeStep)  {
    // Proportional Controller
    Mdouble pController = pGain * stressError;

    //Integral Controller
    iController += iGain * timeStep * stressError;

    // Controller Command
    return pController + iController;
}


