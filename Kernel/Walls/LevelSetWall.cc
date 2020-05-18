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

#include <limits>

#include "LevelSetWall.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"
#include "InfiniteWall.h"

// Constructor; currently only allows predefined shapes
LevelSetWall::LevelSetWall(Shape s, double radius, ParticleSpecies* sp) : radius_(radius)
{
    setSpecies(sp);
    if (s == Shape::Sphere)
    {
        setShapeSphere();
        createVTKSphere();
    }
    else if (s == Shape::Cube)
    {
        setShapeCube();
        createVTKCube();
    }
    else if (s == Shape::Diamond)
    {
        setShapeDiamond();
        createVTK();
    }
    else if (s == Shape::Cylinder)
    {
        setShapeCylinder();
        createVTK();
    }
    else if (s == Shape::FourSided)
    {
        setShapeFourSided();
        createVTK();
    }
    else
    {
        logger(ERROR, "Shape unknown");
    }
}

LevelSetWall::~LevelSetWall()
{
    logger(DEBUG, "LevelSetWall::~LevelSetWall finished");
}

/*!
 * Wall copy method. It calls the copy constructor of this Wall, useful for polymorphism
 */
LevelSetWall* LevelSetWall::copy() const
{
    return new LevelSetWall(*this);
}

///*!
// * \param[in] otherPosition The position to which the distance must be computed to.
// * \return The distance of the wall to the particle.
// */
//Mdouble LevelSetWall::getDistance(Vec3D otherPosition) const
//{
//    return getOrientation().getDistance(otherPosition,getPosition());
//}

/*!
 * \param[in] p BaseParticle for which the distance to the wall must be computed.
 * \param[out] distance Distance between the particle and the wall.
 * \param[out] normal_return The normal of this wall, will only be set if there is a collision.
 * \return A boolean value for whether or not there is a collision.
 * \details First the distance is checked. If there is no collision, this
 * function will return false and give the distance. If there is a collision, the
 * function will return true and give the distance and the normal vector of this wall.
 * Since this function should be called before calculating any 
 * Particle-Wall interactions, it can also be used to set the normal vector in 
 * case of curved walls.
 */
bool LevelSetWall::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    //transform coordinates into position-orientation frame
    Vec3D position = p.getPosition() - getPosition();
    getOrientation().rotateBack(position);
    ///\todo do this for all walls
    if (getDistanceAndNormalLabCoordinates(position, p.getWallInteractionRadius(this), distance, normal_return))
    {
        getOrientation().rotate(normal_return);
        return true;
    }
    return false;
    //return getDistanceAndNormal(p.getPosition(), distance, normal_return, p.getRadius(), p.getInteractionRadius());
}

/*!
 * \param[in] is The input stream from which the LevelSetWall is read. Only needed for backward compatibility.
 */
void LevelSetWall::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy;
    ///\todo
}

/*!
 * \return The string "LevelSetWall", which is the name of this class.
 */
std::string LevelSetWall::getName() const
{
    return "LevelSetWall";
}

void LevelSetWall::writeVTK(VTKContainer& vtk) const
{
    vtk.triangleStrips = vtkLabFrame_.triangleStrips;
    vtk.points = vtkLabFrame_.points;
    for (Vec3D& p : vtk.points)
    {
        getOrientation().rotate(p);
        p += getPosition();
    }
}

void LevelSetWall::createVTKCube()
{
    Mdouble r = radius_ * (N - 1) / N;
    vtkLabFrame_.points.resize(8);
    vtkLabFrame_.points[0] = r * Vec3D(-1, -1, -1);
    vtkLabFrame_.points[1] = r * Vec3D(-1, -1, 1);
    vtkLabFrame_.points[2] = r * Vec3D(-1, 1, -1);
    vtkLabFrame_.points[3] = r * Vec3D(-1, 1, 1);
    vtkLabFrame_.points[4] = r * Vec3D(1, -1, -1);
    vtkLabFrame_.points[5] = r * Vec3D(1, -1, 1);
    vtkLabFrame_.points[6] = r * Vec3D(1, 1, -1);
    vtkLabFrame_.points[7] = r * Vec3D(1, 1, 1);
    vtkLabFrame_.triangleStrips.resize(6);
    vtkLabFrame_.triangleStrips[0] = {0, 1, 2, 3};
    vtkLabFrame_.triangleStrips[1] = {4, 5, 6, 7};
    vtkLabFrame_.triangleStrips[2] = {0, 1, 4, 5};
    vtkLabFrame_.triangleStrips[3] = {2, 3, 6, 7};
    vtkLabFrame_.triangleStrips[4] = {0, 2, 4, 6};
    vtkLabFrame_.triangleStrips[5] = {1, 3, 5, 7};
}

void LevelSetWall::createVTKDiamond()
{
    vtkLabFrame_.points.resize(6);
    vtkLabFrame_.points[0] = radius_ * Vec3D(-1, 0, 0);
    vtkLabFrame_.points[1] = radius_ * Vec3D(1, 0, 0);
    vtkLabFrame_.points[2] = radius_ * Vec3D(0, -1, 0);
    vtkLabFrame_.points[3] = radius_ * Vec3D(0, 1, 0);
    vtkLabFrame_.points[4] = radius_ * Vec3D(0, 0, -1);
    vtkLabFrame_.points[5] = radius_ * Vec3D(0, 0, 1);
    vtkLabFrame_.triangleStrips.resize(8);
    vtkLabFrame_.triangleStrips[0] = {0, 3, 4};
    vtkLabFrame_.triangleStrips[1] = {0, 3, 5};
    vtkLabFrame_.triangleStrips[2] = {0, 2, 4};
    vtkLabFrame_.triangleStrips[3] = {0, 2, 5};
    vtkLabFrame_.triangleStrips[4] = {1, 3, 4};
    vtkLabFrame_.triangleStrips[5] = {1, 3, 5};
    vtkLabFrame_.triangleStrips[6] = {1, 2, 4};
    vtkLabFrame_.triangleStrips[7] = {1, 2, 5};
}

void LevelSetWall::createVTK()
{
    Mdouble distance;
    Vec3D normal;
    InfiniteWall wall;

    Mdouble interactionRadius = radius_/(2*N);
    for (int i=-N; i<N; ++i) {
        for (int j=-N; j<N; ++j) {
            for (int k=-N; k<N; ++k) {
                Vec3D min = interactionRadius*Vec3D(2*i,2*j,2*k);
                Vec3D max = min + interactionRadius*Vec3D(2,2,2);
                Vec3D position = min + interactionRadius*Vec3D(1,1,1);
                if (getDistanceAndNormalLabCoordinates(position, 2*interactionRadius, distance, normal)) {
                    wall.set(normal,position+distance*normal);
                    std::vector<Vec3D> points;
                    wall.createVTK(points,min,max);
                    if (!points.empty()) {
                        wall.addToVTK(points, vtkLabFrame_);
                    }
                }
            }
        }
    }
}

void LevelSetWall::createVTKSphere()
{
    vtkLabFrame_.points.clear();
    vtkLabFrame_.triangleStrips.clear();
    
    const unsigned nr = 40;
    const unsigned nz = 40;
    
    std::array<Mdouble, nr> s, c;
    for (unsigned i = 0; i < nr; ++i)
    {
        s[i] = sin(2.0 * constants::pi * i / nr);
        c[i] = cos(2.0 * constants::pi * i / nr);
    }
    std::array<Mdouble, nz> h, r;
    for (unsigned j = 0; j < nz; ++j)
    {
        r[j] = radius_ * sin(2.0 * constants::pi * j / nz);
        h[j] = radius_ * cos(2.0 * constants::pi * j / nz);
    }
    
    for (unsigned j = 0; j < nz; ++j)
    {
        for (unsigned i = 0; i < nr; ++i)
        {
            vtkLabFrame_.points.emplace_back(r[j] * s[i], r[j] * c[i], h[j]);
        }
    }
    
    vtkLabFrame_.triangleStrips.reserve(nz - 1);
    for (unsigned j = 0; j < nz - 1; ++j)
    {
        std::vector<double> cell;
        cell.reserve(2 * nr + 2);
        for (unsigned i = 0; i < nr; ++i)
        {
            cell.push_back(i + j * nr);
            cell.push_back(i + (j + 1) * nr);
        }
        cell.push_back(j * nr);
        cell.push_back((j + 1) * nr);
        vtkLabFrame_.triangleStrips.push_back(cell);
    }
}

void LevelSetWall::writeToFile(int n, double extraRadius) const
{
    std::stringstream s;
    s << "Values of level set function\n";
    Vec3D normal;
    Mdouble distance;
    double meshSize = (radius_ + extraRadius) / n;
    for (int k = -n; k <= n; ++k)
    {
        for (int j = -n; j <= n; ++j)
        {
            for (int i = -n; i <= n; ++i)
            {
                Vec3D position = meshSize * Vec3D(i, j, k);
                getDistanceAndNormalLabCoordinates(position, radius_ + extraRadius, distance, normal);
                s << position << ' '
                  << distance << ' '
                  << normal << '\n';
            }
        }
    }
    helpers::writeToFile("LevelSet.txt", s.str());
    helpers::writeToFile("LevelSet.m", "close all\n"
                                       "data=importdata('LevelSet.txt');\n"
                                       "N = nthroot(size(data.data,1),3);\n"
                                       "x = reshape(data.data(:,1),N,N,N);\n"
                                       "y = reshape(data.data(:,2),N,N,N);\n"
                                       "z = reshape(data.data(:,3),N,N,N);\n"
                                       "d = reshape(data.data(:,4),N,N,N);\n"
                                       "nx= reshape(data.data(:,5),N,N,N);\n"
                                       "ny= reshape(data.data(:,6),N,N,N);\n"
                                       "nz= reshape(data.data(:,7),N,N,N);\n"
                                       "mid = (N-1)/2;\n"
                                       "surf(x(:,:,mid),y(:,:,mid),d(:,:,mid),'FaceColor','none')\n"
                                       "hold on\n"
                                       "quiver(x(:,:,mid),y(:,:,mid),nx(:,:,mid),ny(:,:,mid))\n"
                                       "legend('distance','normal');\n"
                                       "if (min(min(d(:,:,mid)))<0 && max(max(d(:,:,mid)))>0)\n"
                                       "    contourf(x(:,:,mid),y(:,:,mid),d(:,:,mid),[0 0])\n"
                                       "    legend('distance','normal','distance>=0'); set(legend,'Location','best');\n"
                                       "end\n"
                                       "hold off\n"
                                       "xlabel('x'); ylabel('y'); title('cross-section at z=0');\n"
                                       "view(0,0); axis equal;");
    logger(INFO, "Run LevelSet.m to view level set");
}

bool LevelSetWall::getDistanceAndNormalLabCoordinates(Vec3D position, Mdouble interactionRadius, Mdouble& distance,
                                                      Vec3D& normal) const
{
    //scale values [-radius,radius] to [0,2*N]
    Vec3D positionScaled = position / radius_ * static_cast<double>(N) + Vec3D(1, 1, 1) * static_cast<double>(N);
    // index of level set value used for interpolation (the discrete level-set value left of the interpolation point)
    int i = std::max(std::min((int) positionScaled.X, 2 * N - 1), 0);
    int j = std::max(std::min((int) positionScaled.Y, 2 * N - 1), 0);
    int k = std::max(std::min((int) positionScaled.Z, 2 * N - 1), 0);
    // scaled difference between position and discrete level-set value, in [0,1) for interpolation
    double x = positionScaled.X - i;
    double y = positionScaled.Y - j;
    double z = positionScaled.Z - k;
    // trilinear interpolation
    distance = levelSet_[i][j][k] * (1 - x) * (1 - y) * (1 - z)
               + levelSet_[i + 1][j][k] * x * (1 - y) * (1 - z)
               + levelSet_[i][j + 1][k] * (1 - x) * y * (1 - z)
               + levelSet_[i + 1][j + 1][k] * x * y * (1 - z)
               + levelSet_[i][j][k + 1] * (1 - x) * (1 - y) * z
               + levelSet_[i + 1][j][k + 1] * x * (1 - y) * z
               + levelSet_[i][j + 1][k + 1] * (1 - x) * y * z
               + levelSet_[i + 1][j + 1][k + 1] * x * y * z;
    if (distance < interactionRadius)
    {
        normal.X = -levelSet_[i][j][k] * (1 - y) * (1 - z)
                   + levelSet_[i + 1][j][k] * (1 - y) * (1 - z)
                   - levelSet_[i][j + 1][k] * y * (1 - z)
                   + levelSet_[i + 1][j + 1][k] * y * (1 - z)
                   - levelSet_[i][j][k + 1] * (1 - y) * z
                   + levelSet_[i + 1][j][k + 1] * (1 - y) * z
                   - levelSet_[i][j + 1][k + 1] * y * z
                   + levelSet_[i + 1][j + 1][k + 1] * y * z;
        normal.Y = -levelSet_[i][j][k] * (1 - x) * (1 - z)
                   - levelSet_[i + 1][j][k] * x * (1 - z)
                   + levelSet_[i][j + 1][k] * (1 - x) * (1 - z)
                   + levelSet_[i + 1][j + 1][k] * x * (1 - z)
                   - levelSet_[i][j][k + 1] * (1 - x) * z
                   - levelSet_[i + 1][j][k + 1] * x * z
                   + levelSet_[i][j + 1][k + 1] * (1 - x) * z
                   + levelSet_[i + 1][j + 1][k + 1] * x * z;
        normal.Z = -levelSet_[i][j][k] * (1 - x) * (1 - y)
                   - levelSet_[i + 1][j][k] * x * (1 - y)
                   - levelSet_[i][j + 1][k] * (1 - x) * y
                   - levelSet_[i + 1][j + 1][k] * x * y
                   + levelSet_[i][j][k + 1] * (1 - x) * (1 - y)
                   + levelSet_[i + 1][j][k + 1] * x * (1 - y)
                   + levelSet_[i][j + 1][k + 1] * (1 - x) * y
                   + levelSet_[i + 1][j + 1][k + 1] * x * y;
        Mdouble length = normal.getLength();
        if (length != 0)
            normal /= -length;
        else
            normal = {1,0,0};
        return true;
    }
    return false;
}

void LevelSetWall::setShapeSphere()
{
    logger(INFO, "Creating a sphere");
    for (int k = -N; k <= N; ++k)
    {
        for (int j = -N; j <= N; ++j)
        {
            for (int i = -N; i <= N; ++i)
            {
                levelSet_[i + N][j + N][k + N] = radius_ * sqrt(i * i + j * j + k * k) / N - radius_;
            }
        }
    }
}

void LevelSetWall::setShapeCylinder()
{
    logger(INFO, "Creating a cylinder");
    for (int k = -N; k <= N; ++k)
    {
        for (int j = -N; j <= N; ++j)
        {
            for (int i = -N; i <= N; ++i)
            {
                levelSet_[i + N][j + N][k + N] = radius_*std::fmax(sqrt(i * i + j * j) / N - 1,fabs(k)/N-1);
            }
        }
    }
}

void LevelSetWall::setShapeCube()
{
    logger(INFO, "Creating a cube");
    radius_ *= N / (N - 1);
    Mdouble h = radius_ / N;
    for (int k = -N; k <= N; ++k)
    {
        for (int j = -N; j <= N; ++j)
        {
            for (int i = -N; i <= N; ++i)
            {
                std::array<int, 3> sort;
                sort[0] = abs(i);
                sort[1] = abs(j);
                sort[2] = abs(k);
                std::sort(sort.begin(), sort.end());
                if (sort[1] < N)
                {
                    // face
                    levelSet_[i + N][j + N][k + N] = h * (sort[2] - (N - 1));
                }
                else if (sort[0] < N)
                {
                    // edge
                    levelSet_[i + N][j + N][k + N] = h * 2;
                }
                else
                {
                    //vertex
                    levelSet_[i + N][j + N][k + N] = h * 3;
                }
                //edge and vertex distances longer than real, otherwise trilinear extrapolation is problematic.
            }
        }
    }
}

void LevelSetWall::setShapeDiamond()
{
    logger(INFO, "Creating a diamond square");
    for (int k = -N; k <= N; ++k)
    {
        for (int j = -N; j <= N; ++j)
        {
            for (int i = -N; i <= N; ++i)
            {
                levelSet_[i + N][j + N][k + N] = radius_ / sqrt(3) / N * ((double) (abs(i) + abs(j) + abs(k)) - N);
            }
        }
    }
}

void LevelSetWall::setShapeFourSided()
{
    //open filestream
    std::string fileName = "fourSided.txt";
    std::ifstream is(fileName.c_str(), std::ios::in);
    if (is.fail())
    {
        helpers::writeToFile("fourSided.m",
                "%% parameters\n"
                "n=10; % 2n+1 is the number of cells in the level set per dimension\n"
                "w=1; % [-w,w] is the width of the shape in z-direction\n"
                "% define a few points on the xy-cross section of the shape that will be interpolated\n"
                "% we define a quarter of the cross section only; thus the first point has differ \n"
                "% from the last point by a quarter rotation around the origin\n"
                "x0=[1 1.8 1.8 1.4 1];\n"
                "y0=[-1 0 1 1.4 1];\n"
                "\n"
                "%% create curve, and plot\n"
                "% define a parameter t0, such that the parametrisation is x0(t0), y0(t0)\n"
                "t0=1:length(x0);\n"
                "% define a more fine-grained parameter t1, and interpolate, obtaining a \n"
                "% finer parametrisation x1(t1), y1(t1)\n"
                "t1= linspace(1,length(x0));\n"
                "x1=interpn(t0,x0,t1,'spline');\n"
                "y1=interpn(t0,y0,t1,'spline');\n"
                "% complete the shape by repeating it four times, each rotated by 90 deg\n"
                "% the points (x,y) are on the shape \n"
                "x=[x1 -y1 -x1 y1];\n"
                "y=[y1 x1 -y1 -x1];\n"
                "% plot the resulting shape\n"
                "figure(1)\n"
                "plot(x,y,'-',x0,y0,'o')\n"
                "axis equal\n"
                "\n"
                "%% create level set, and plot\n"
                "% compute the size of the domain needed to contain the shape ...\n"
                "radius = max([max(abs(x)),max(abs(y)),abs(w)]);\n"
                "% .. then divide by the number of cells to get the cell size\n"
                "radius = n/floor(n/radius);\n"
                "% create xLS,yLS,zLS values of the level set function dist(xLS,yLS,zLS)\n"
                "LS = linspace(-radius,radius,2*n+1);\n"
                "[xLS,yLS,zLS]=meshgrid(LS);\n"
                "% define the level set function by computing the distance of (xLS,yLS,zLS)\n"
                "% to the nearest point on the shape (x,y,+-w)\n"
                "dist = zeros(size(xLS));\n"
                "for i=1:length(xLS(:))\n"
                "   dist(i) = sqrt(min((x-xLS(i)).^2+(y-yLS(i)).^2));\n"
                "   % make it a signed distance\n"
                "   if inpolygon(xLS(i),yLS(i),x,y)\n"
                "       dist(i) = -dist(i);\n"
                "   end\n"
                "   % compute the signed distance as the max of signed distances in xy and z\n"
                "   dist(i) = max(dist(i),abs(zLS(i))-w);\n"
                "end\n"
                "% rescale the distance\n"
                "dist = dist/radius;\n"
                "\n"
                "% plot the i-th cross-section in xy (i=n+1 is in the center)\n"
                "figure(2)\n"
                "i=n+4\n"
                "mesh(xLS(:,:,i),yLS(:,:,i),dist(:,:,i),'FaceColor','none')\n"
                "hold on\n"
                "contourf(xLS(:,:,i),yLS(:,:,i),dist(:,:,i),[0 0])\n"
                "plot(x0,y0,'ro')\n"
                "hold off\n"
                "view(0,90)\n"
                "axis equal\n"
                "title(['z=' num2str(unique(zLS(:,:,i)))])\n"
                "xlabel('x'), ylabel('y'), zlabel('dist')\n"
                "\n"
                "% plot the i-th cross-section in yz (i=n+1 is in the center)\n"
                "figure(3)\n"
                "mesh(squeeze(yLS(:,i,:)),squeeze(zLS(:,i,:)),squeeze(dist(:,i,:)),'FaceColor','none')\n"
                "hold on\n"
                "contourf(squeeze(yLS(:,i,:)),squeeze(zLS(:,i,:)),squeeze(dist(:,i,:)),[0 0])\n"
                "hold off\n"
                "view(90,0)\n"
                "axis equal\n"
                "title(['x=' num2str(unique(xLS(:,i,:)))])\n"
                "xlabel('y'), ylabel('z'), zlabel('dist')\n"
                "\n"
                "% plot the i-th cross-section in xz (i=n+1 is in the center)\n"
                "figure(4)\n"
                "mesh(squeeze(xLS(i,:,:)),squeeze(zLS(i,:,:)),squeeze(dist(i,:,:)),'FaceColor','none')\n"
                "hold on\n"
                "contourf(squeeze(xLS(i,:,:)),squeeze(zLS(i,:,:)),squeeze(dist(i,:,:)),[0 0])\n"
                "hold off\n"
                "view(90,0)\n"
                "axis equal\n"
                "title(['y=' num2str(unique(yLS(i,:,:)))])\n"
                "xlabel('x'), ylabel('z'), zlabel('dist')\n"
                "\n"
                "%% write n and dist values to file (read into Mercury)\n"
                "fid = fopen('fourSided.txt','w');\n"
                "fprintf(fid,'%d\\n',n);\n"
                "fprintf(fid,'%f ',dist);\n"
                "fclose(fid);");
        logger(ERROR, "readFromFile: file % could not be opened; make sure you run fourSided.m first.", fileName);
    }
    
    int n;
    is >> n;
    if (n != N)
    {
        logger(ERROR, "readFromFile: level set size does not match %", n);
    }
    
    Mdouble unitDistance;
    for (int k = -N; k <= N; ++k)
    {
        for (int j = -N; j <= N; ++j)
        {
            for (int i = -N; i <= N; ++i)
            {
                is >> unitDistance;
                levelSet_[i + N][j + N][k + N] = radius_ * unitDistance;
            }
        }
    }
}
