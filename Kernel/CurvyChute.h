#ifndef CURVYCHUTE_H
#define CURVYCHUTE_H
#include "Mercury3D.h"
#include "Particles/BaseParticle.h"

#include<functional>

/*!
 * \class CurvyChute
 * \brief Creates chutes defined by curvilinear coordinates. Inherits from Mercury3D. 
 * \details Like the Chute class, we get three new features compared to the
 * parent Mercury3D class: 
 * <ul>
 *   <li> setChuteAngle </li>
 * </ul>
 *
 */

class CurvyChute : public Mercury3D
{
    public:
        CurvyChute();

        ~CurvyChute();

        Mdouble getChuteAngle() const;

        Mdouble getChuteAngleDegrees() const;

        void setChuteAngleAndMagnitudeOfGravity(Mdouble chuteAngle, Mdouble gravity);

        void setSurface( std::function<Vec3D(Mdouble ccR, Mdouble ccS)> surface,
                         std::function<Mdouble(Mdouble ccR, Mdouble ccS)> areaElement,
                         Mdouble ccRMin, Mdouble ccRMax, Mdouble ccdR, 
                         Mdouble ccSMin, Mdouble ccSMax, Mdouble ccdS );

        void setBasalThickness(Mdouble basalThickness);

        void setBasalDensity(Mdouble basalDensity);

        void setBasalDisorder(Mdouble basalDisorder);

        void setBasalSizeDispersity(Mdouble basalSizeDispersity);

        /*!
         * \brief The basal prototype should already have a species assigned to
         * it. 
         */
        void setBasalPrototype(BaseParticle* basalPrototype);

        /*!
         * \brief When this is called, a copy of basalPrototype->getSpecies()
         * will be added to the speciesHandler, and assigned to basalSpecies_.
         * (The original species might already be there, which is fine.) 
         */
        void createBottom();

        /*!
         * \brief When this is called, all particles in the particleHandler that
         * have the species basalSpecies_ will be deleted, and createBottom()
         * will be called again.
         */
        void recreateBottom();

    private:
        /*! 
         * \brief Chute angle. 
         */
        Mdouble chuteAngle_;

        Mdouble magnitudeOfGravity_;

        /*! 
         * \brief The shape of the chute, defined by a function of two
         * curvilinear coordinates. The 'r' direction is the streamwise
         * coordinate while the 's' coordinate is the cross-stream coordinate.
         *
         * \details JMFT: By convention, gravity will point in the 'x'
         * and 'z' directions, and the chute should be considered flat.
         */
        std::function<Vec3D(Mdouble ccR, Mdouble ccS)> surface_;
        
        /*!
         * \brief The area element of the surface
         */
        std::function<Mdouble(Mdouble ccR, Mdouble ccS)>  areaElement_;

        Mdouble ccRMin_, ccRMax_, ccdR_, ccSMin_, ccSMax_, ccdS_;

#if 0
        /*!
         * \brief Scale factor (length of d(surface)/dr), as a function of r and s
         */
        std::function<Mdouble(Mdouble ccR, Mdouble ccS)> scaleFactorR;

        /*!
         * \brief Scale factor (length of d(surface)/ds), as a function of r and s
         */
        std::function<Mdouble(Mdouble ccR, Mdouble ccS)> scaleFactorS;
#endif


        /*! 
         * \brief The thickness of the base, in terms of number of layers
         */
        Mdouble basalThickness_;

        /*!
         * \brief The amount of disorder in the base, from 0 (ordered) to 1
         * (completely disordered)
         */
        Mdouble basalDisorder_;

        /*!
         * \brief The area density of the basal particles.
         */
        Mdouble basalDensity_;

        /*!
         * \brief The amount of size dispersity amongst the basal particles.
         */
        Mdouble basalSizeDispersity_;

        /*!
         * \brief Pointer to a prototypical basal particle (not in a
         * particleHandler, but should have a species set)
         */
        BaseParticle* basalPrototype_;

        /*!
         * \brief Pointer to the species of the prototypical basal particle.
         */
        ParticleSpecies* basalSpecies_; 

       /*!
        * \brief Determines whether the chute should have periodic (TRUE) or
        * solid (FALSE) walls in the s direction. 
        * \todo JMFT: How to reconcile this with possibly curvilinear
        * coordinates, curvilinear edges?
        */
        // bool isChutePeriodic_;
};

#endif
