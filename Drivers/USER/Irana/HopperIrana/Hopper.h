//
// Created by irana on 9/13/17.
//

#ifndef MERCURY_HOPPER
#define MERCURY_HOPPER

#include "WallHandler.h"
#include "BaseObject.h"

/*!
 * Class that creates a hopper with the given input parameters. After setting the parameters, makeHopper() can be called
 * to make the 2 or 4 IntersectionOfWalls that together form the hopper, and add them to the wallHandler.
 */
class Hopper
{
public:
    
    Hopper();
    
    ///\todo implement write, read, copy
    Hopper* copy() const;
    void read(std::istream& is);
    void write(std::ostream& os) const;
    
    std::string getName() const
    {
        return "Hopper";
    }
    
    Mdouble getHopperLength() const;
    
    void setHopperLength(Mdouble hopperLength_);
    
    Mdouble getHopperHeight() const;
    
    void setHopperHeight(Mdouble hopperHeight_);
    
    Mdouble getHopperAngle() const;
    
    void setHopperAngle(Mdouble hopperAngle_);
    
    Mdouble getHopperExitLength() const;
    
    void setHopperExitLength(Mdouble hopperExitLength_);
    
    Mdouble getHopperShift() const;
    
    void setHopperShift(Mdouble hopperShift_);
    
    Mdouble getHopperLift() const;
    
    void setHopperLift(Mdouble hopperLift_);
    
    unsigned int getHopperDimension() const;
    
    void setHopperDimension(unsigned int hopperDimension_);
    
    Mdouble getHopperLowestPoint() const;
    
    void setHopperLowestPoint(Mdouble hopperLowestPoint_);
    
    void makeHopper(WallHandler& wallHandler);

private:
    
    /*!
     * \brief Dimension of the hopper in vertical direction
     */
    Mdouble hopperLength_;
    
    /*!
     * \brief Dimension of the hopper in horizontal direction
     */
    Mdouble hopperHeight_;
    
    /*!
     * \brief Angle between the two pieces of the hopper walls
     */
    Mdouble hopperAngle_;
    
    /*!
     * \brief Dimension of the hopper exit in vertical direction
     */
    Mdouble hopperExitLength_;
    
    /*!
     * \brief The x position where the Chute starts (defined as the beginning of the hopper)
     */
    Mdouble hopperShift_;
    
    /*!
     * \brief This is the vertical distance the chute is lifted above the plane.
     */
    Mdouble hopperLift_;
    
    /*!
     * \brief This is the dimension of the hopper, by default it is one dimensional and hence does not have side wall
     */
    unsigned int hopperDimension_;
    
    /*!
     * \brief The z coordinate of the right C' point (when the left C point is in the origin)
     */
    Mdouble hopperLowestPoint_;
};


#endif //MERCURY_HOPPER
