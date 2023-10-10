

#include "../../../../oomph-lib/src/meshes/simple_cubic_mesh.h"
#include "Elements/AndersonJackson.h"

template <class ELEMENT>
class UnderResolvedMesh : public oomph::Problem
{
public:
    UnderResolvedMesh() {};
    ~UnderResolvedMesh() {};

    UnderResolvedMesh(const double& xMin, const double& xMax,
                      const double& yMin, const double& yMax,
                      const double& zMin, const double& zMax,
                      const unsigned int& nx, const unsigned int& ny, const unsigned int& nz,
                      oomph::TimeStepper* time_stepper_pt = &oomph::Mesh::Default_TimeStepper)
    {

    };

    /// Doc the solution
    virtual void doc_solution(oomph::DocInfo& doc_info_) {}
    /// Doc the voidage
    virtual void doc_voidage(oomph::DocInfo& doc_info_) {}
    /// Doc the elements
    virtual void doc_element(oomph::DocInfo& doc_info_) {}

    oomph::DocInfo doc_info;

};

