#include "../../../../oomph-lib/src/meshes/simple_cubic_mesh.h"

template<class ELEMENT>
class UnderResolvedMesh : public oomph::RefineableSimpleCubicMesh<ELEMENT>
{
public:
    UnderResolvedMesh()= default;

    UnderResolvedMesh(const double& xMin, const double& xMax,
                      const double& yMin, const double& yMax,
                      const double& zMin, const double& zMax,
                      const unsigned int& nx, const unsigned int& ny, const unsigned int& nz,
                      oomph::TimeStepper* time_stepper_pt = &oomph::Mesh::Default_TimeStepper)
    {
        oomph::Problem::mesh_pt() = new oomph::RefineableSimpleCubicMesh<ELEMENT>(xMin, xMax, yMin, yMax, zMin, zMax, nx, ny, nz, time_stepper_pt);
    };

};

