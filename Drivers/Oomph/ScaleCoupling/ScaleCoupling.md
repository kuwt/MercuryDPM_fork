# Volume Coupling in MercuryDPM

This directory contains three codes:
- ParticleBeam: A DPM simulation using MercuryDPM,
- SolidBeam: A FEM simulation using oomph-lib,
- ScaleCoupledBeam: A coupled simulation using Moomph.

All three codes solve the same problem: a solid beam that is excited by a fixed velocity on the left end and a free boundary at the right end. The FEM elements are the same size as the DEM particles.

The DEM and FEM simulations agree well, see
https://youtu.be/hbL-c1EHK6g.


# How ScaleCoupling works

```plantuml
class ScaleCoupledBeam {}

ScaleCoupledBeam --|> ScaleCoupledSolid : ELEMENT = RefineableQDPVDElement<3, 2>

class ScaleCoupledSolid<ELEMENT>

ScaleCoupledSolid --|> ScaleCoupling : SELEMENT = ScaleCoupledElement<ELEMENT>

class ScaleCoupling<Mercury3D,SolidProblem<SELEMENT>> {
std::function couplingWeight_
std::vector<CoupledElement> coupledElements_
std::vector<CoupledParticle> coupledParticles_
double penalty_
}

ScaleCoupling --|> BaseCoupling

class CoupledElement {
ELEMENT* el_pt;
std::vector<BaseParticle*> coupledParticles
}

ScaleCoupling --> CoupledElement : uses

class CoupledParticle {
    Vec3D couplingForce
    CoupledElement* coupledElement_pt
    Vector<double> s
}

ScaleCoupling --> CoupledParticle : uses

class ScaleCoupledElement<ELEMENT> {
    Vector<Vector<double>> coupling_residual
    Vector<Vector<double>> coupling_jacobian
    Vector<double> coupling_weight
}

ScaleCoupling -> ScaleCoupledElement : uses

BaseCoupling --|> Mercury3D

class Mercury3D {}

class BaseCoupling<Mercury3D,SolidProblem<SELEMENT>> {}

BaseCoupling --|> SolidProblem

class SolidProblem<SELEMENT> {
    std::string name_
    double elasticModulus_
    double poissonRatio_
    double density_
    double gravity_
    function_pt body_force_fct
    ConstitutiveLaw* constitutive_law_pt
    SolidMesh* Solid_mesh_pt
    SolidMesh* Traction_mesh_pt
    std::function isPinned_
}
```

#How SCoupling works

```plantuml
class SCoupledDriver {
  CoupledBeam() : sets up mercury, oomph, calls solveSurfaceCoupling()
}

class SCoupledSolidProblem<ELEMENT> {
}

class SCoupling<M,O> {
 void solveSurfaceCoupling() : calls solveSurfaceCoupling(dtFEM/dtDEM)
 void solveSurfaceCoupling(int nStep) : calls getSCoupledElements(), \n\tcreateDPMWallsFromFiniteElems(), computeOneTimeStepForSCoupling()
 getSCoupledElements()
 createDPMWallsFromFiniteElems()
 computeOneTimeStepForSCoupling(nStep) : calls updateDPMWallsFromFiniteElems(), \n\tsolveMercury(), updateTractionOnFiniteElems(), solveOomph()
 updateDPMWallsFromFiniteElems() 
 updateTractionOnFiniteElems()
}

class BaseCoupling<M,O> {
 void solveMercury(nSteps) : calls Mercury3D::computeOneTimeStep()
 void solveOomph() : calls actionsBeforeOomphTimeStep() and unsteady_newton_solve()
}

class SolidProblem<ELEMENT> {
 std::string name_
 double elasticModulus_
 double poissonRatio_
 double density_
 double gravity_
 function body_force_fct
 ConstitutiveLaw* constitutive_law_pt
 SolidMesh* Solid_mesh_pt
 SolidMesh* Traction_mesh_pt
 std::function isPinned_
 virtual void actionsBeforeOomphTimeStep()
 virtual void actionsBeforeSolve()
 virtual void actionsAfterSolve()
}

SCoupledDriver --|> SCoupledSolidProblem : ELEMENT = RefineableQDPVDElement<3, 2>
SCoupledSolidProblem --|> SCoupling : M = Mercury3D, \nO = SolidProblem<SCoupledElement<ELEMENT>>
SCoupling --|> BaseCoupling
BaseCoupling --|> SolidProblem
```

# How SCoupling should work

```plantuml
class ScaleCoupledBeam {
  CoupledBeam() : calls solveVolumeCoupling()
}

class ScaleCoupledSolidProblem<ELEMENT> {
}

class ScaleCoupling<M,O> {
 Vector<DPMScaleCoupledElement> listOfDPMScaleCoupledElement;
 void solveVolumeCoupling() : \n\tcalls solveVolumeCoupling(dtFEM/dtDEM)
 void solveVolumeCoupling(int nStep) : \n\tcalls initialise, computeTimeStep
 void initialise() : calls M::initialiseSolve, O::prepareForSolve
 void computeTimeStep(nStep) : \n\tcalls updateOomph,O::solveUnsteady, \n\tupdateMercuryDPM, M::computeTimeStep
 void updateOomph()
 void updateMercuryDPM()
}

class SolidBeam {
SolidBeam() : calls \n\tprepareForSolve,\n\tsolveUnsteady
}

class SolidProblem<ELEMENT> {
 std::string name_
 double elasticModulus_
 double poissonRatio_
 double density_
 double gravity_
 function body_force_fct
 ConstitutiveLaw* constitutive_law_pt
 SolidMesh* Solid_mesh_pt
 SolidMesh* Traction_mesh_pt
 std::function isPinned_
 void prepareForSolve()
 void solveUnsteady(timeMax, dt, saveCount)
 virtual void actionsBeforeOomphTimeStep()
 virtual void actionsBeforeSolve()
 virtual void actionsAfterSolve()
}

class ParticleBeam {
ParticleBeam() : \n\tcalls solve()
}

class Mercury3D {
 void solve() : calls \n\tinitialiseSolve, \n\tcomputeTimeStep, \n\tfinaliseSolve
 void initialiseSolve()
 void computeTimeStep()
 void finaliseSolve()
}

class ScaleCoupledElement<ELEMENT> {
 Vector<Vector<double> > nodal_coupling_residual;
 Vector<Vector<double> > nodal_coupling_jacobian;
 Vector<Vector<double> > couplingStiffness;
}

class CoupledSolidNode
{
double coupling_weight;
Vector<double> coupling_force;
}

class DPMScaleCoupledElement
{
ELEMENT* bulk_elem_pt;
Vector<BaseParticle*> listOfCoupledParticles;
Vector<BaseParticle*> listOfCoupledParticlesExt;
Vector<Vector<double>> listOfParticleCentersLoc;
Vector<Vector<double>> listOfShapesAtParticleCenters;
Vector<Vector<double>> couplingMatrix;
}
    
ParticleBeam --|> Mercury3D 
SolidBeam --|> SolidProblem : ELEMENT = RefineableQDPVDElement<3, 2>
ScaleCoupledBeam --|> ScaleCoupledSolidProblem : ELEMENT = RefineableQDPVDElement<3, 2>
ScaleCoupledSolidProblem --|> ScaleCoupling : M = Mercury3D, \nO = SolidProblem<ScaleCoupledElement<ELEMENT>>
ScaleCoupling --|> SolidProblem
ScaleCoupling --> ScaleCoupledElement : uses
ScaleCoupledElement --> CoupledSolidNode : uses
ScaleCoupling --> DPMScaleCoupledElement : uses
ScaleCoupling --|> Mercury3D
```
Todo's:
- remove Mercury3D as template parameters
