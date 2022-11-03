//
// Created by reza on 11-01-21.
//

#include "PController.h"
PController::PController(Mdouble pGain)
: pGain(pGain)
        {}

Mdouble PController::apply (Mdouble stressError, Mdouble timeStep) {
// Proportional Controller
    Mdouble pController = pGain * stressError;

// Controller Command
    return pController;
}