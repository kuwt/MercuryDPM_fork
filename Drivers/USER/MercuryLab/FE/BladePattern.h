#ifndef MERCURY_BLADEPATTERN_H
#define MERCURY_BLADEPATTERN_H

#include "Logger.h"

enum class BladePattern : unsigned {
    LongScrew, //OldScrew long
    ShortScrew, //OldScrew short
    FullMixing, //3
    FullMixingWithForward, //3F
    FullMixingWithForwardWithPins, //3F*
    FullMixingShifted, //8
    FullAlt, //6
    FullAltShifted, //11
};

/*!
 * write to file
 */
std::ostream& operator<<(std::ostream& os, BladePattern type)
{
    switch (type) {
        case BladePattern::LongScrew: os << "LongScrew"; break;
        case BladePattern::ShortScrew: os << "ShortScrew"; break;
        case BladePattern::FullMixing: os << "FullMixing"; break;
        case BladePattern::FullMixingWithForward: os << "FullMixingWithForward"; break;
        case BladePattern::FullMixingWithForwardWithPins: os << "FullMixingWithForwardWithPins"; break;
        case BladePattern::FullMixingShifted: os << "FullMixingShifted"; break;
        case BladePattern::FullAlt: os << "FullAlt"; break;
        case BladePattern::FullAltShifted: os << "FullAltShifted"; break;
        default: logger(ERROR,"BladePattern % not found", static_cast<unsigned>(type));
    }
    return os;
}


//set blade pattern
std::string getBladePatternString(BladePattern pattern)
{
    //note, pattern starts at outflow, ends at inflow, unlike in the Excel sheet provided by FE
    switch (pattern)
    {
        case BladePattern::LongScrew: //original screw pattern with forward blades in the middle
            return " f   b   b   b   f   f   f   f  "
                   "   f   f   f   f   f   b   b   f"
                   " f   b   b   b   f   f   f   f  "
                   "   f   f   f   f   f   b   b   f";
        case BladePattern::ShortScrew:
            return "   f   f   f  "
                   " f   b   b   f"
                   "   f   f   f  "
                   " f   b   b   f";
        case BladePattern::FullMixing: //simplified blade configuration: forward-backward blade combinations
            return " f   b   b   b   b   b   f"
                   "   f   f   f   f   f   f  "
                   " f   b   b   b   b   b   f"
                   "   f   f   f   f   f   f  ";
        case BladePattern::FullMixingWithForward: //added forward blades in the middle
            return " f   b   b   f   b   b   f"
                   "   f   f   f   f   f   f  "
                   " f   b   b   f   b   b   f"
                   "   f   f   f   f   f   f  ";
        case BladePattern::FullMixingWithForwardWithPins: //added pins for deagglomeration
            return " f   b   b   f   b   b   f"
                   "   f   f   fp pf  pf  pf  "
                   " f   b   b   f   b   b   f"
                   "   f   f   f   f   f   f  ";
        case BladePattern::FullMixingShifted: //two center rows are shifted by one
            return " f   b   b   b   b   b   f"
                   "  f   f   f   f   f   f   "
                   "f   b   b   b   b   b   f "
                   "   f   f   f   f   f   f  ";
        case BladePattern::FullAlt: //forward-forward-backward-backward
            return " f   b   b   b   b   b   f"
                   "   f   f   f   f   f   f  "
                   " f   f   f   f   f   f   f"
                   "   f   b   b   b   b   f  ";
        case BladePattern::FullAltShifted: //two center rows are shifted by one
            return " f   b   b   b   b   b   f"
                   "  f   f   f   f   f   f   "
                   "f   f   f   f   f   f   f "
                   "   f   b   b   b   b   f  ";
        default:
            logger(ERROR,"This blade pattern does not exist:", static_cast<unsigned>(pattern));
    }
}

Mdouble getBladeNumber(BladePattern pattern) {
    std::string patternString = getBladePatternString(pattern);
    return patternString.length()/8;
}


#endif //MERCURY_BLADEPATTERN_H
