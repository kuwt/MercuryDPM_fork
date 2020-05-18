//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Particles/SuperQuadricParticle.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Mercury3D.h"

///Tests whether the shape-function of a superquadric and its gradient and hessian are computed correctly
///For now, only test spheres and ellipsoids.
class ShapeGradientHessianTester : public Mercury3D
{
public:
    void test()
    {
        testSphere();
        testEllipsoid();
        testRoundedBeam();
        testCushion();
    
        logger(INFO, "All tests pass.");
    }
    
    void testSphere()
    {
        SuperQuadricParticle p;
        auto species = new LinearViscoelasticSpecies();
        speciesHandler.addObject(species);
        p.setSpecies(species);
        //sphere of radius 1
        p.setAxesAndExponents(1, 1, 1, 1, 1);
        Vec3D point = Vec3D(1, 0, 0);
        logger.assert_always(mathsFunc::isEqual(p.computeShape(point), 0, 1e-5),
                             "shape function of sphere should be 0 at (1,0,0)");
        logger.assert_always(mathsFunc::isEqual(p.computeShapeGradientLabFixed(point), Vec3D(2, 0, 0), 1e-5),
                             "gradient of shapefunction sphere wrong");
        Matrix3D hessian;
        hessian.XX = 2;
        hessian.YY = 2;
        hessian.ZZ = 2;
        logger.assert_always(mathsFunc::isEqual(p.computeHessianLabFixed(point), hessian, 1e-5),
                             "hessian of sphere wrong");
    
    
        point = sqrt(1. / 3) * Vec3D(1, 1, 1);
        logger.assert_always(mathsFunc::isEqual(p.computeShape(point), 0, 1e-5),
                             "shape function of sphere should be 0 at 1/sqrt(3) * (1,1,1)");
        logger.assert_always(mathsFunc::isEqual(p.computeShapeGradientLabFixed(point),
                                                2.0 / std::sqrt(3) * Vec3D(1, 1, 1), 1e-5),
                             "gradient of shapefunction sphere wrong");
        logger.assert_always(mathsFunc::isEqual(p.computeHessianLabFixed(point), hessian, 1e-5),
                             "hessian of sphere wrong");
    
        point = Vec3D(0, 0, 0);
        logger.assert_always(p.computeShape(point) < 0,
                             "Origin should be inside sphere, so negative shape value");
    
        point = Vec3D(2, 0, 0);
        logger.assert_always(p.computeShape(point) > 0, "Shape sphere should be positive outside sphere");
    }
    
    void testEllipsoid()
    {
        SuperQuadricParticle p;
        auto species = new LinearViscoelasticSpecies();
        speciesHandler.addObject(species);
        p.setSpecies(species);
        p.setAxesAndExponents(3, 2, 1, 1, 1);
        p.setOrientationViaNormal({0.0, 1.0, 0.0});
        Vec3D point = Vec3D(2, 0, 0);
        logger.assert_always(mathsFunc::isEqual(p.computeShape(point), 0, 1e-5),
                             "shape function of rotated ellipsoid should be 0 at (2,0,0)");
        Vec3D gradientBodyFixed = p.computeShapeGradientLabFixed(Vec3D(2, 0, 0));
        p.getOrientation().rotateBack(gradientBodyFixed);
        logger.assert_always(mathsFunc::isEqual(gradientBodyFixed, Vec3D(0, -1, 0), 1e-5),
                             "lab-fixed gradient of shapefunction ellipsoid wrong, is % should be (0,-1,0)",
                             gradientBodyFixed, point);
        Matrix3D hessian;
        hessian.XX = 2. / 9;
        hessian.YY = 1. / 2;
        hessian.ZZ = 2;
        SmallMatrix<3, 3> hessianLabFixed = p.computeHessianLabFixed(point);
        SmallMatrix<3, 3> rotationMatrix;
        p.getOrientation().getRotationMatrix(rotationMatrix);
        SmallMatrix<3, 3> hessianBodyFixed = rotationMatrix.transpose() * (hessianLabFixed) * rotationMatrix;
        logger.assert_always(mathsFunc::isEqual(hessianBodyFixed, hessian, 1e-5), "hessian of ellipsoid wrong");
    
        point = sqrt(1. / 3) * Vec3D(2, 3, 1);
        logger.assert_always(mathsFunc::isEqual(p.computeShape(point), 0, 1e-5),
                             "shape function of rotated ellipsoid should be 0 at std::sqrt(1./3) * (2,3,1)");
        gradientBodyFixed = p.computeShapeGradientLabFixed(sqrt(1. / 3) * Vec3D(2, 3, 1));
        p.getOrientation().rotateBack(gradientBodyFixed);
        logger.assert_always(mathsFunc::isEqual(gradientBodyFixed, sqrt(1. / 3) * Vec3D(2. / 3, -1, 2), 1e-5),
                             "body-fixed gradient of shapefunction ellipsoid wrong, is % should be "
                                     "std::sqrt(1./3) (2./3, -1, 2)", gradientBodyFixed);
    
        hessianLabFixed = p.computeHessianLabFixed(point);
        p.getOrientation().getRotationMatrix(rotationMatrix);
        hessianBodyFixed = rotationMatrix.transpose() * (hessianLabFixed) * rotationMatrix;
        logger.assert_always(mathsFunc::isEqual(hessianBodyFixed, hessian, 1e-5), "hessian of ellipsoid wrong");
    
        logger.assert_always(p.computeShape(Vec3D(0, 0, 0)) < 0,
                             "Origin should be inside ellipsoid, so negative shape value");
        logger.assert_always(p.computeShape(Vec3D(3, 0, 0)) > 0, "Shape ellipsoid should be positive outside");
    }
    
    void testRoundedBeam()
    {
        SuperQuadricParticle p;
        auto species = new LinearViscoelasticSpecies();
        speciesHandler.addObject(species);
        p.setSpecies(species);
        p.setAxesAndExponents(3, 2, 1, 0.5, 0.5);
        Vec3D point = Vec3D(3, 2, 1);
        logger.assert_always(mathsFunc::isEqual(p.computeShape(point), 2, 1e-5),
                             "shape function of rounded beam should be 2 at (3,2,1)");
        Vec3D gradientBodyFixed = p.computeShapeGradientLabFixed(point);
        p.getOrientation().rotateBack(gradientBodyFixed);
        logger.assert_always(mathsFunc::isEqual(gradientBodyFixed, Vec3D(4.0/3.0, 2.0, 4.0), 1e-5),
                             "lab-fixed gradient of shapefunction rounded beam wrong", gradientBodyFixed, point);
        Matrix3D hessian;
        hessian.XX = 4. / 3;
        hessian.YY = 3;
        hessian.ZZ = 12.0;
        SmallMatrix<3, 3> hessianLabFixed = p.computeHessianLabFixed(point);
        SmallMatrix<3, 3> rotationMatrix;
        p.getOrientation().getRotationMatrix(rotationMatrix);
        SmallMatrix<3, 3> hessianBodyFixed = rotationMatrix.transpose() * (hessianLabFixed) * rotationMatrix;
        logger.assert_always(mathsFunc::isEqual(hessianBodyFixed, hessian, 1e-5), "hessian of rounded beam wrong");
    }
    
    void testCushion()
    {
        SuperQuadricParticle p;
        auto species = new LinearViscoelasticSpecies();
        speciesHandler.addObject(species);
        p.setSpecies(species);
        p.setAxesAndExponents(3, 2, 1, 0.5, 1);
        Vec3D point = Vec3D(-1, 1, 1);
        logger.assert_always(mathsFunc::isEqual(p.computeShape(point), (1./9 + 1./4.0)* (1./9 + 1./4.0), 1e-5),
                             "shape function of cushion wrong  at (-1,1,1)");
        Vec3D gradientBodyFixed = p.computeShapeGradientLabFixed(point);
        p.getOrientation().rotateBack(gradientBodyFixed);
        logger.assert_always(mathsFunc::isEqual(gradientBodyFixed, Vec3D(-13.0/81.0, 13.0/36.0, 4.0), 1e-5),
                             "lab-fixed gradient of shapefunction cushion wrong", gradientBodyFixed, point);
        Matrix3D hessian;
        hessian.XX = 7. / 27;
        hessian.YY = 31. / 36;
        hessian.ZZ = 12.0;
        hessian.XY = -2.0/9.0;
        hessian.YX = hessian.XY;
        SmallMatrix<3, 3> hessianLabFixed = p.computeHessianLabFixed(point);
        SmallMatrix<3, 3> rotationMatrix;
        p.getOrientation().getRotationMatrix(rotationMatrix);
        SmallMatrix<3, 3> hessianBodyFixed = rotationMatrix.transpose() * (hessianLabFixed) * rotationMatrix;
        logger.assert_always(mathsFunc::isEqual(hessianBodyFixed, hessian, 1e-5), "hessian of cushion wrong");
    }
};

int main()
{
    ShapeGradientHessianTester test;
    test.test();
    return 0;
}
