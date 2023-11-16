#include "CurvyChute.h"

CurvyChute::CurvyChute()
{
    chuteAngle_ = 0;
    surface_ = [] (Mdouble ccR, Mdouble ccS) { return Vec3D(ccR, ccS, 0); } ;
    ccRMin_ = 0;
    ccRMax_ = 1;
    ccdR_ = 0.05;
    ccSMin_ = 0;
    ccSMax_ = 1;
    ccdS_ = 0.05;

    basalThickness_ = 1;
    basalDensity_ = 1;
    basalDisorder_ = 1;
    basalSizeDispersity_ = 0;

    basalPrototype_ = NULL;
    basalSpecies_ = NULL;
}

CurvyChute::~CurvyChute() {
    /*
    if (basalPrototype_ != NULL)
        delete basalPrototype_;

    if (basalSpecies_ != NULL)
        delete basalSpecies_;
    */
}

Mdouble CurvyChute::getChuteAngle() const
{
    return chuteAngle_;
}

Mdouble CurvyChute::getChuteAngleDegrees() const
{
    return getChuteAngle() * 180. / constants::pi;
}  

void CurvyChute::setChuteAngleAndMagnitudeOfGravity(Mdouble chuteAngle, Mdouble magnitudeOfGravity)
{
    chuteAngle_ = chuteAngle;
    magnitudeOfGravity_ = magnitudeOfGravity;
    setGravity(
            Vec3D(
                magnitudeOfGravity_ * sin(chuteAngle_),
                0,
               -magnitudeOfGravity_ * cos(chuteAngle_)
                ));
}

void CurvyChute::setSurface( std::function<Vec3D(Mdouble ccR, Mdouble ccS)> surface,
        std::function<Mdouble(Mdouble ccR, Mdouble ccS)> areaElement,
        Mdouble ccRMin, Mdouble ccRMax, Mdouble ccdR, 
        Mdouble ccSMin, Mdouble ccSMax, Mdouble ccdS )
{
    if (ccdR <= 0 || ccdS <= 0)
        logger(ERROR, "[CurvyChute::setSurface] ccdR % and ccdS % should both be positive.", 
                ccdR, ccdS);

    surface_     = surface;
    areaElement_ = areaElement;
    ccRMin_      = ccRMin;
    ccRMax_      = ccRMax;
    ccdR_        = ccdR;
    ccSMin_      = ccSMin;
    ccSMax_      = ccSMax;
    ccdS_        = ccdS;
}

void CurvyChute::setBasalThickness(Mdouble basalThickness)
{

    basalThickness_ = basalThickness;
}

void CurvyChute::setBasalDensity(Mdouble basalDensity)
{
    if (basalDensity < 0)
        logger(ERROR, "[CurvyChute::setBasalDensity(Mdouble)] basal density % is negative", basalDensity);
    basalDensity_ = basalDensity;
}

void CurvyChute::setBasalDisorder(Mdouble basalDisorder)
{
    if (basalDisorder < 0)
        logger(ERROR, "[CurvyChute::setBasalDisorder(Mdouble)] basal disorder % is negative", basalDisorder);
    basalDisorder_ = basalDisorder;
}

void CurvyChute::setBasalSizeDispersity(Mdouble basalSizeDispersity)
{
    if (basalSizeDispersity >= 0 && basalSizeDispersity < 1) 
        basalSizeDispersity_ = basalSizeDispersity;
    else
        logger(ERROR, "[CurvyChute::setBasalSizeDispersity(Mdouble)] size dispersity % is not in [0,1)",
                basalSizeDispersity);
}

void CurvyChute::setBasalPrototype(BaseParticle* basalPrototype)
{
    basalPrototype_ = basalPrototype->copy();
}

void CurvyChute::createBottom()
{
    /* Make a new species for the base. (Doesn't matter if the original species
     * was already in the speciesHandler, make a new one.) */
    basalSpecies_ = speciesHandler.copyAndAddObject(
            basalPrototype_->getSpecies()->copy());
    /* Update the species of the basalPrototype_ */
    basalPrototype_->setSpecies( basalSpecies_ );  

    if (basalThickness_ < 0 || basalDensity_ == 0)
    {
        logger(WARN, "[CurvyChute::createBottom()] basalThickness %, basalDensity % imply that no chute is necessary!",
                basalThickness_, basalDensity_);
        return;
    }

    if (ccRMin_ >= ccRMax_ || ccSMin_ >= ccSMax_)
        logger(ERROR, "[CurvyChute::createBottom()] The bounds of the coordinates are not in order.");

    for (Mdouble r = ccRMin_ + ccdR_ * basalDisorder_ * random.getRandomNumber(-1, 1); 
            r <= ccRMax_; r += ccdR_)
    {
        for (Mdouble s = ccSMin_ + ccdS_ * basalDisorder_ * random.getRandomNumber(-1, 1); 
                s <= ccSMax_; s += ccdS_)
        {
            Mdouble area = areaElement_(r, s) * ccdR_ * ccdS_;
            Mdouble areaOfOneParticle = constants::pi * pow(basalPrototype_->getRadius(), 2);
            unsigned int numberToInsert = random.getPoissonVariate( basalDensity_ * area / areaOfOneParticle ); 
            // unsigned int numberToInsert = floor( basalDensity_ * area / areaOfOneParticle ); 

            for (int j = 0; j < numberToInsert; j++)
            {
                auto particleToAdd = basalPrototype_->copy();
                particleToAdd->setRadius( 
                        basalPrototype_->getRadius()  
                            * (1 + basalSizeDispersity_ * random.getRandomNumber(-1,1)) );
                Mdouble rpos = r + random.getRandomNumber(-.5, .5) * basalDisorder_ * ccdR_;
                Mdouble spos = s + random.getRandomNumber(-.5, .5) * basalDisorder_ * ccdS_;
                particleToAdd->setPosition(surface_(rpos, spos) + Vec3D(0, 0, random.getRandomNumber(-.5, .5) * basalThickness_));
                particleToAdd->setVelocity(Vec3D(0,0,0));
                particleToAdd->fixParticle();
                particleHandler.addObject(particleToAdd);
                // delete particleToAdd;
            }
        }
    }
}

/*!
 * \details Do not run this in the middle of a simulation: This does not check
 * for interactions with existing particles, so you could create a basal
 * particle that intersects and interacts strongly with one. 
 */
void CurvyChute::recreateBottom()
{
    for (auto p = particleHandler.begin(); p != particleHandler.end(); ++p)
    {
        if ((*p)->getSpecies() == basalSpecies_)
        {
            particleHandler.removeObject((*p)->getId());
            p--;
        }
    }
    speciesHandler.removeObject(basalSpecies_->getId()); 
    createBottom();
}
