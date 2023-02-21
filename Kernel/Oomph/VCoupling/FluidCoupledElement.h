//
// Created by mitchel on 11/7/22.
//

#ifndef FLUIDCOUPLEDELEMENT_H
#define FLUIDCOUPLEDELEMENT_H

namespace oomph
{

//=================start_wrapper==================================
/// Wrapper class for fluid elements to be coupled with discrete
/// solid particles in a overlapping volume
//================================================================
    template<class ELEMENT>
    class FluidCoupledElement : public virtual ELEMENT
    {

    public:
    
        /// Constructor: Call constructor of underlying element
        FluidCoupledElement()
        {};
    
        /// Destructor (empty)
        ~FluidCoupledElement()
        {};
    
        void setFluidVolumeFraction(double fvf)
        {
            if (fvf != 1.0) {logger(DEBUG,"Set fluid volume fraction to %",fvf);}
            fluidVolumeFraction = fvf;
        }
        double getFluidVolumeFaction() {return fluidVolumeFraction;}
        
    private:
        double fluidVolumeFraction; //Default to 1.0 in case of a non-UnderResolved method
    
    };
}


#endif //FLUIDCOUPLEDELEMENT_H
