/* Parabolic chute
 * Based on the Coil class */

#ifndef PARABOLA_H
#define PARABOLA_H

#include "BaseWall.h"
#include "Math/Vector.h"

class ParabolaChute : public BaseWall {
    public:
        /*!
         * \brief Default constructor, sets a chute with default parameters.
         */
        ParabolaChute();
        /*!
         * \brief Copy constructor 
         */
        ParabolaChute(const ParabolaChute& other);
        /*!
         * \brief Constructor in which all parameters are set.
         */
        ParabolaChute(Mdouble length, Mdouble widthscale);
        ~ParabolaChute();
        void set(Mdouble length, Mdouble widthscale);
        ParabolaChute* copy() const override;
        bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const override;
        std::vector<BaseInteraction*> getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) override;
        void read(std::istream& is) override;
        void write(std::ostream& os) const override;
        std::string getName() const override;
        
    private:
        Mdouble l_;
        Mdouble ws_;
};

#endif
