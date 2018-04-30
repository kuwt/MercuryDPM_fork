/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// This plugin has been developed and contributed by Sven Buijssen, TU
// Dortmund, Germany.
// Thanks to Bryn Lloyd (blloyd@vision.ee.ethz.ch) at ETH Zuerich for
// developing and sharing vtkTensorGlyphFilter, the ancestor of this
// filter. That filter's output (i.e. spheres) can be mimicked by setting both
// ThetaRoundness and PhiRoundness to 1.0.
// Thanks to Gordon Kindlmann for pointing out at VisSym04 that superquadrics
// can be put to a good use to visualize tensors and for pointing out to me
// an insufficient initial implementation of this filter.

#include "vtkSuperquadricTensorGlyphFilter.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTensorGlyphSameEigensystem.h"
#include "vtkSuperquadricSource.h"
#include "vtkCellArray.h"
#include "vtkAppendPolyData.h"
#include "vtkMath.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"


vtkStandardNewMacro(vtkSuperquadricTensorGlyphFilter);

vtkSuperquadricTensorGlyphFilter::vtkSuperquadricTensorGlyphFilter()
{
  this->SetNumberOfInputPorts(1);
  this->ThetaResolution = 16;
  this->PhiResolution = 16;
  this->ThetaRoundness = 0.3;
  this->PhiRoundness = 0.3;
  this->Gamma = 1.5;
  this->ScaleFactor = 1.0;
  this->ExtractEigenvalues = 0;

  // by default, process active point scalars
  //this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
  //                             vtkDataSetAttributes::SCALARS);

  // by default, process active point scalars
  //this->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
  //                            vtkDataSetAttributes::VECTORS);

  // by default, process active point tensors
  //this->SetInputArrayToProcess(2, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
  //                            vtkDataSetAttributes::TENSORS);

  // by default, process active point species type
  this->SetInputArrayToProcess(1, 0, 0, 0,"speciesType");

  // by default, process active point phiAndTheta
  this->SetInputArrayToProcess(2, 0, 0, 0,"phiAndTheta");

  // by default, process active point axes scales
  this->SetInputArrayToProcess(3, 0, 0, 0,"axesScales");

  // by default, process active point orientation tensor
  this->SetInputArrayToProcess(4, 0, 0, 0,"orientationTensor");

}

vtkSuperquadricTensorGlyphFilter::~vtkSuperquadricTensorGlyphFilter()
{
}

int vtkSuperquadricTensorGlyphFilter::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}


int vtkSuperquadricTensorGlyphFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPointData *inputPD;
  vtkPointData *outputPD;
  vtkSuperquadricSource *superquadric;
  vtkAppendPolyData* append;
  vtkIdType numPtsIn, numGlyphPtsTotal, ptIncr, inPtId, i;
  int j, axis, abort = 0;
  vtkTensorGlyphSameEigensystem** tensorsArray = NULL;
  vtkIdType tensorsArrayCount = 0;
  vtkTensorGlyphSameEigensystem* tensors = NULL;

  vtkDataArray *inSpeciesType = NULL;
  vtkDataArray *inPhiAndTheta = NULL;  
  vtkDataArray *inAxesScales = NULL;  
  vtkDataArray *inOrientationTensor = NULL;

  double speciesType[1];
  double phiAndTheta[2];
  double axesScales[3];
  double tensor[9];

  double *m[3], w[3], *v[3];
  double m0[3], m1[3], m2[3];
  double v0[3], v1[3], v2[3];
  vtkDoubleArray *lambda;
  vtkIdTypeArray *numGlyphPts;
  double sum, cl, cp, maxScale, theta, phi; // cs
  vtkIdType id[1];
  int showWarning = 0;

  // set up working matrices
  m[0] = m0; m[1] = m1; m[2] = m2;
  v[0] = v0; v[1] = v1; v[2] = v2;

  inputPD = input->GetPointData();
  if (!inputPD)
    {
    vtkErrorMacro(<<"No data to glyph!");
    return 1;
    }

  numPtsIn = input->GetNumberOfPoints();
  if (numPtsIn < 1)
    {
    vtkErrorMacro(<<"No points to glyph!");
    return 1;
    }

  // Store eigenvalues for later useas data array
  if ( this->ColorGlyphs &&
       this->ColorMode == COLOR_BY_EIGENVALUES)
    {
    lambda = vtkDoubleArray::New();
    lambda->SetNumberOfComponents(3);
    lambda->SetNumberOfTuples(numPtsIn);
    lambda->SetName("Eigenvalues");
    }

  inSpeciesType = this->GetInputArrayToProcess(1, inputVector);
  inPhiAndTheta = this->GetInputArrayToProcess(2, inputVector);
  inAxesScales = this->GetInputArrayToProcess(3, inputVector);
  inOrientationTensor = this->GetInputArrayToProcess(4, inputVector);

  if (!inOrientationTensor)
    {
    vtkErrorMacro(<<"No data to glyph!");
    return 1;
    }

  append = vtkAppendPolyData::New();

  // Not every input point necessarily has valid tensor data and is
  // hence glyphed; additionally, not every glyph necessarily has the
  // same number of points.
  numGlyphPts = vtkIdTypeArray::New();
  numGlyphPts->SetNumberOfValues(numPtsIn);
  tensorsArray = new vtkTensorGlyphSameEigensystem*[numPtsIn];
  for (inPtId = 0; inPtId < numPtsIn && !abort; inPtId++)
    {
    if (! (inPtId % 500))
      {
      // this->UpdateProgress(static_cast<double>(inPtId)/numPtsIn);
      // abort = this->GetAbortExecute();
      }

    if (!abort)
      {     
      
      inSpeciesType->GetTuple(inPtId, speciesType);           
      inPhiAndTheta->GetTuple(inPtId, phiAndTheta);
      inAxesScales->GetTuple(inPtId, axesScales); 
      inOrientationTensor->GetTuple(inPtId, tensor);    

      // Init eigenvalues
      w[0] = w[1] = w[2] = 0.0;

        // Copy point inPtId and its point data to a single point input object
        vtkPolyData* singlePoint = vtkPolyData::New();
        vtkPoints* newPt = vtkPoints::New();
        vtkCellArray* newVertex = vtkCellArray::New();
        newPt->SetNumberOfPoints(1);
        newPt->SetPoint(0,input->GetPoint(inPtId));
        singlePoint->SetPoints(newPt);
        newPt->Delete();
        singlePoint->GetPointData()->CopyAllOn();
        singlePoint->GetPointData()->CopyAllocate(inputPD, 1);
        singlePoint->GetPointData()->CopyData(inputPD, inPtId, 0);

        // Create vertex for this point
        id[0] = inPtId;
        newVertex->InsertNextCell(1,id);
        singlePoint->SetVerts(newVertex);
        newVertex->Delete();

        // set the active tensor for next operation
        singlePoint->GetPointData()->SetActiveTensors( inOrientationTensor->GetName() );
        if (inSpeciesType)
          {
          singlePoint->GetPointData()->SetActiveScalars( inSpeciesType->GetName() );
          }

        tensors = vtkTensorGlyphSameEigensystem::New();
        tensors->SetInputData(singlePoint);
        tensors->SetExtractEigenvalues(ExtractEigenvalues);
        tensors->SetColorGlyphs(this->GetColorGlyphs());
        tensors->SetColorMode(this->GetColorMode());
        tensors->SetScaleFactor(ScaleFactor);
        tensors->SetClampScaling(ClampScaling);
        tensors->SetMaxScaleFactor(MaxScaleFactor);
        tensors->SetThreeGlyphs(0);
        singlePoint->Delete();

        superquadric = vtkSuperquadricSource::New();
        superquadric->SetThetaResolution(this->ThetaResolution);
        superquadric->SetPhiResolution(this->PhiResolution);
        superquadric->ToroidalOff();
        axis = 2;  // superquadrics with rotational symmetriy around the z axis

        // Set roundness per tensor, if requested
        this->PhiRoundness = phiAndTheta[0];
        this->ThetaRoundness = phiAndTheta[1];

        if (!this->FixedThetaPhiRoundness)
          {
          theta = phiAndTheta[0];
          phi = phiAndTheta[1];

          w[0] = axesScales[0];
          w[1] = axesScales[1];
          w[2] = axesScales[2];

          v[0][0] = tensor[0]; v[1][0] = tensor[3]; v[2][0] = tensor[6];
          v[0][1] = tensor[1]; v[1][1] = tensor[4]; v[2][1] = tensor[7];
          v[0][2] = tensor[2]; v[1][2] = tensor[5]; v[2][2] = tensor[8];

          // avoid calculating them again in vtkTensorGlyph:
          tensors->SetEigenvectorsEigenValues(v, w);
      }

        superquadric->SetAxisOfSymmetry(axis);
        superquadric->SetThetaRoundness(theta);
        superquadric->SetPhiRoundness(phi);
        superquadric->Update();
//        tensors->SetSourceConnection(superquadric->GetOutputPort());
        tensors->SetSourceData(superquadric->GetOutput());
        tensors->SetScaleFactor(this->ScaleFactor);
        tensors->Update();
        numGlyphPts->SetValue(inPtId, tensors->GetOutput()->GetNumberOfPoints());
        superquadric->Delete();

        // Add tensor glyph filter's output to output
        // fast variant
        tensorsArray[tensorsArrayCount++] = tensors;
        // slow variant
//        append->AddInputData( vtkPolyData::SafeDownCast( tensors->GetOutput()) );
//        append->Update();
//        tensors->Delete();
       
      }

    if ( this->ColorGlyphs &&
         this->ColorMode == COLOR_BY_EIGENVALUES)
      {
      maxScale = 1.0;
      if ( this->ClampScaling )
        {
        for (maxScale=0.0, i=0; i<3; i++)
          {
          if ( maxScale < fabs(w[i]) )
            {
            maxScale = fabs(w[i]);
            }
          }
        if ( maxScale > this->MaxScaleFactor )
          {
          maxScale = this->MaxScaleFactor / maxScale;
          }
        }
      for (i=0; i<3; i++)
        {
        lambda->SetComponent(inPtId, i, w[i] * maxScale);
        }
      }
    } // end for loop

  append->SetUserManagedInputs(1);
  append->SetNumberOfInputs(tensorsArrayCount);
  vtkIdType count;
  for (count = 0; count < tensorsArrayCount && !abort; count++)
    {
    tensors = tensorsArray[count];
    if (tensors != NULL)
      {
      append->SetInputConnectionByNumber(count, tensors->GetOutputPort());
      }
    }
  append->Update();
  for (count = 0; count < tensorsArrayCount && !abort; count++)
    {
    tensors = tensorsArray[count];
    if(tensors != NULL)
      {
      tensors->Delete();
      tensorsArray[count] = NULL;
      tensors = NULL;
      }
    }
  delete [] tensorsArray;
  tensorsArray = NULL;

  if (append)
    output->ShallowCopy(append->GetOutput());

  if (showWarning)
    {
    vtkWarningMacro(<< "Tensors with non-positive eigenvalues have been omitted. This "
                    << "filter assumes three real, positive eigenvalues.");
    }

  outputPD = output->GetPointData();

  // In constrast to vtkTensorGlyph vtkTensorGlyphSameEigensystem names the
  // data fields it creates such that ParaView will color them. No need
  // to compensate here.

  numGlyphPtsTotal = 0;
  for (inPtId = 0; inPtId < numPtsIn; inPtId++)
    {
      numGlyphPtsTotal += numGlyphPts->GetValue(inPtId);
    }

  // Add eigenvalues as data array
  if ( this->ColorGlyphs &&
       this->ColorMode == COLOR_BY_EIGENVALUES)
    {
    vtkDoubleArray *eigenvalues = vtkDoubleArray::New();
    eigenvalues->SetNumberOfComponents(3);
    eigenvalues->SetNumberOfTuples(numGlyphPtsTotal);
    eigenvalues->SetName(lambda->GetName());
    ptIncr = 0;
    for (inPtId=0; inPtId < numPtsIn; inPtId++)
      {
      for (i=0; i < numGlyphPts->GetValue(inPtId); i++)
        {
        eigenvalues->SetTuple(ptIncr+i, inPtId, lambda);
        }
      ptIncr += numGlyphPts->GetValue(inPtId);
      }
    int idx = outputPD->AddArray(eigenvalues);
    outputPD->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    lambda->Delete();
    eigenvalues->Delete();

    // Remove data array named MaxEigenvalue automatically created
    // by vtkTensorGlyphSameEigensystem as that data array is comprised in
    // the vector data array named Eigenvalues provided by this class.
    if (outputPD->GetArray("MaxEigenvalue"))
      {
      outputPD->RemoveArray("MaxEigenvalue");
      }
    }
  numGlyphPts->Delete();
  append->Delete();

  return 1;
}


void vtkSuperquadricTensorGlyphFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
