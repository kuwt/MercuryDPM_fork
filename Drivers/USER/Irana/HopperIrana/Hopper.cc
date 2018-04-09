#include <Walls/IntersectionOfWalls.h>
#include "Hopper.h"
#include "WallHandler.h"
#include "DPMBase.h"

Hopper::Hopper()
{
    hopperLength_ = 0;
    hopperHeight_ = 0;
    hopperAngle_ = 0;
    hopperExitLength_ = 0;
    hopperShift_ = 0;
    hopperLift_ = 0;
    hopperLowestPoint_ = 0;
    hopperDimension_ = 1;
}

Mdouble Hopper::getHopperLength() const
{
    return hopperLength_;
}

void Hopper::setHopperLength(Mdouble hopperLength)
{
    hopperLength_ = hopperLength;
}

Mdouble Hopper::getHopperHeight() const
{
    return hopperHeight_;
}

void Hopper::setHopperHeight(Mdouble hopperHeight)
{
    hopperHeight_ = hopperHeight;
}

Mdouble Hopper::getHopperAngle() const
{
    return hopperAngle_;
}

void Hopper::setHopperAngle(Mdouble hopperAngle)
{
    hopperAngle_ = hopperAngle;
}

Mdouble Hopper::getHopperExitLength() const
{
    return hopperExitLength_;
}

void Hopper::setHopperExitLength(Mdouble hopperExitLength)
{
    hopperExitLength_ = hopperExitLength;
}

Mdouble Hopper::getHopperShift() const
{
    return hopperShift_;
}

void Hopper::setHopperShift(Mdouble hopperShift)
{
    hopperShift_ = hopperShift;
}

Mdouble Hopper::getHopperLift() const
{
    return hopperLift_;
}

void Hopper::setHopperLift(Mdouble hopperLift)
{
    hopperLift_ = hopperLift;
}

unsigned int Hopper::getHopperDimension() const
{
    return hopperDimension_;
}

void Hopper::setHopperDimension(unsigned int hopperDimension)
{
    hopperDimension_ = hopperDimension;
}

Mdouble Hopper::getHopperLowestPoint() const
{
    return hopperLowestPoint_;
}

void Hopper::setHopperLowestPoint(Mdouble hopperLowestPoint)
{
    hopperLowestPoint_ = hopperLowestPoint;
}

void Hopper::makeHopper(WallHandler& wallHandler)
{
    //hopper walls
    //to create the finite hopper walls, we take vector between two wall points in xz-plane, then rotate clockwise and make unit length
    // A\       /A'
    //   \     /       A,B,C denote three points on the left and right hopper walls which are used to construct the hopper
    //    \   /        shift denotes the space by which the chute has to be shifted to the right such that the hopper is in the domain
    //   B|   |B'
    //    |   |
    //   C|   |
    //        |C'
    
    Vec3D gravityDirection = wallHandler.getDPMBase()->getGravity();
    gravityDirection /= gravityDirection.getLength();
    const Mdouble s = gravityDirection.X;
    const Mdouble c = -gravityDirection.Z;
    
    // "0.5*(hopperLength_+hopperExitLength_) / tan(hopperAngle_)" is the minimum heigth of the hopper, to make sure things should flow down and not to the sides.
    // hopperHeight_ is now an input variable
    // hopperHeight_ = hopperLowestPoint_ + 1.1 * 0.5*(hopperLength_+hopperExitLength_) / tan(hopperAngle_);
    
    const Mdouble HopperCornerHeight = hopperHeight_ - 0.5 * (hopperLength_ - hopperExitLength_) / std::tan(hopperAngle_);
    // first we create the LEFT hopper wall
    
    // coordinates of A,B,C in (vertical parallel to flow, vertical normal to flow, horizontal) direction
    Vec3D A = Vec3D(0.0, 0.0, hopperHeight_);
    Vec3D B = Vec3D(0.5 * (hopperLength_ - hopperExitLength_), 0.0, HopperCornerHeight);
    Vec3D C = Vec3D(0.5 * (hopperLength_ - hopperExitLength_), 0.0, 0.0);
    
    // now rotate the coordinates of A,B,C to be in (x,y,z) direction, such that C-B is parallel to gravity
    A = Vec3D(c * A.X - s * A.Z, 0.0, s * A.X + c * A.Z);
    B = Vec3D(c * B.X - s * B.Z, 0.0, s * B.X + c * B.Z);
    C = Vec3D(c * C.X - s * C.Z, 0.0, s * C.X + c * C.Z);
    
    A.X += hopperShift_;
    B.X += hopperShift_;
    C.X += hopperShift_;
    
    //This lifts the hopper a distance above the chute
    A.Z += hopperLift_;
    B.Z += hopperLift_;
    C.Z += hopperLift_;
    
    //create a finite wall from B to A and from C to B on the left hand side
    Vec3D temp = B - A;
    Vec3D normal = Vec3D(temp.Z, 0.0, -temp.X) / std::sqrt(temp.getLengthSquared());
    IntersectionOfWalls intersectionOfWalls;
    intersectionOfWalls.createPrism({A,B,C});
    wallHandler.copyAndAddObject(intersectionOfWalls);
    intersectionOfWalls.clear();
    
    //next, do the same for the right wall
    A = Vec3D(hopperLength_, 0.0, hopperHeight_);
    B = Vec3D((hopperLength_ - hopperExitLength_) / 2 + hopperExitLength_, 0.0, HopperCornerHeight);
    C = Vec3D((hopperLength_ - hopperExitLength_) / 2 + hopperExitLength_, 0.0, hopperLowestPoint_);
    
    //This rotates the right points A', B', C'
    A = Vec3D(c * A.X - s * A.Z , 0.0, s * A.X + c * A.Z);
    B = Vec3D(c * B.X - s * B.Z , 0.0, s * B.X + c * B.Z);
    C = Vec3D(c * C.X - s * C.Z, 0.0, s * C.X + c * C.Z);
    
    A.X += hopperShift_;
    B.X += hopperShift_;
    C.X += hopperShift_;
    
    //This lifts the hopper a distance above the chute
    A.Z += hopperLift_;
    B.Z += hopperLift_;
    C.Z += hopperLift_;
    
    //create a finite wall from B to A and from C to B on the right hand side
    intersectionOfWalls.createPrism({A,B,C});
    wallHandler.copyAndAddObject(intersectionOfWalls);
    intersectionOfWalls.clear();
    
    // if hopperDimension_ == 2, create inclined hopper walls (like in the X-direction) also in the Y-direction.
    // (Else, place vertical (possibly periodic) walls in Y-direction somewhere else in your code. )
    if (hopperDimension_ == 2)
    {
        //coordinates of A',B',C' in (vertical parallel to flow,vertical normal to flow, horizontal) direction
        A = Vec3D(0.0, hopperLength_, hopperHeight_);
        B = Vec3D(0.0, (hopperLength_ - hopperExitLength_) / 2 + hopperExitLength_, HopperCornerHeight);
        C = Vec3D(0.0, (hopperLength_ - hopperExitLength_) / 2 + hopperExitLength_, 0.0);
        
        //now rotate the coordinates of A,B,C to be in (x,y,z) direction
        A = Vec3D(c * A.X - s * A.Z, A.Y, s * A.X + c * A.Z);
        B = Vec3D(c * B.X - s * B.Z, B.Y, s * B.X + c * B.Z);
        C = Vec3D(c * C.X - s * C.Z, C.Y, s * C.X + c * C.Z);
        // the position of A determines shift and zmax
        A.X += hopperShift_;
        B.X += hopperShift_;
        C.X += hopperShift_;
        
        //This lifts the hopper a distance above the chute
        A.Z += hopperLift_;
        B.Z += hopperLift_;
        C.Z += hopperLift_;
        
        //create a finite wall from B' to A' and from C' to B'
        intersectionOfWalls.createPrism({A,B,C});
        wallHandler.copyAndAddObject(intersectionOfWalls);
        intersectionOfWalls.clear();
        
        //Now for the front y-wall
        A = Vec3D(0.0, 0.0, hopperHeight_);
        B = Vec3D(0.0, (hopperLength_ - hopperExitLength_) / 2, HopperCornerHeight);
        C = Vec3D(0.0, (hopperLength_ - hopperExitLength_) / 2, 0.0);
        
        //now rotate the coordinates of A,B,C to be in (x,y,z) direction
        A = Vec3D(c * A.X - s * A.Z, A.Y, s * A.X + c * A.Z);
        B = Vec3D(c * B.X - s * B.Z, B.Y, s * B.X + c * B.Z);
        C = Vec3D(c * C.X - s * C.Z, C.Y, s * C.X + c * C.Z);
        // the position of A determines shift and zmax
        A.X += hopperShift_;
        B.X += hopperShift_;
        C.X += hopperShift_;
        
        //This lifts the hopper a distance above the chute
        A.Z += hopperLift_;
        B.Z += hopperLift_;
        C.Z += hopperLift_;
        
        //create a finite wall from B to A and from C to B
        intersectionOfWalls.createPrism({A,B,C});
        wallHandler.copyAndAddObject(intersectionOfWalls);
        intersectionOfWalls.clear();
    }
}
