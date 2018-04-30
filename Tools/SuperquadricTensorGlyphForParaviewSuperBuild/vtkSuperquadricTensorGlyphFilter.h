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

// .NAME vtkSuperquadricTensorGlyphFilter - scale and orient superquadric glyph according to tensor eigenvalues and eigenvectors

// .SECTION Description
// vtkSuperquadricTensorGlyphFilter is a filter that generates a superquadric
// glyph at every point in the input data set. The glyphs are oriented and
// scaled according to eigenvalues and eigenvectors of "tensor" data of the
// input data set, interpreting the entries of the 3x3 matrix as principal axes
// of the superquadric and their norm as the length of these axes. Set both
// roundness values to 0.0 to get rectangular glyphs, set them to 1.0 to get
// ellipsoidal glyphs, set theta roundness to 1.0 and phi roundness to 0.0 to
// get cylindrical glyphs. Other values lead to superquadric glyphs which are
// in general favorable as they can be distinguished easily for all view
// angles. The Superquadric Tensor Glyph filter operates on any type of data
// set. Its output is polygonal.

// .SECTION Thanks
// This plugin has been developed and contributed by Sven Buijssen, TU
// Dortmund, Germany.
// Thanks to Bryn Lloyd (blloyd@vision.ee.ethz.ch) at ETH Zuerich for
// developing and sharing vtkTensorGlyphFilter, the ancestor of this
// filter. That filter's output (i.e. spheres) can be mimicked by setting both
// ThetaRoundness and PhiRoundness to 1.0.
// Thanks to Gordon Kindlmann for pointing out at VisSym04 that superquadrics
// can be put to a good use to visualize tensors and for pointing out to me
// an insufficient initial implementation of this filter.

// .SECTION See Also
// vtkTensorGlyph

#ifndef vtkSuperquadricTensorGlyphFilter_h
#define vtkSuperquadricTensorGlyphFilter_h

#include "vtkPolyDataAlgorithm.h"

class vtkSuperquadricTensorGlyphFilter : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkSuperquadricTensorGlyphFilter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkSuperquadricTensorGlyphFilter *New();

  // Description:
  // Set the number of points in the longitude direction. Initial value is 16.
  vtkSetMacro(ThetaResolution,int);

  // Description:
  // Set the number of points in the latitude direction. Initial value is 16.
  vtkSetMacro(PhiResolution,int);

  // Description:
  // Set Superquadric east/west roundness.
  // Values range from 0 (rectangular) to 1 (circular) to higher orders.
  // Initial value is 1.0.
  vtkSetMacro(ThetaRoundness,double);

  // Description:
  // Set Superquadric north/south roundness.
  // Values range from 0 (rectangular) to 1 (circular) to higher orders.
  // Initial value is 1.0.
  vtkSetMacro(PhiRoundness,double);

  // Description:
  // Set/Get roundness for superquadric in case settings are chosen per tensor. Initial value is 1.5.
  vtkGetMacro(Gamma, double);
  vtkSetMacro(Gamma, double);

  // Description:
  // If true, then set theta and phi roundness settings for superquadric per tensor. False by default.
  vtkGetMacro(FixedThetaPhiRoundness, int);
  vtkSetMacro(FixedThetaPhiRoundness, int);
  vtkBooleanMacro(FixedThetaPhiRoundness, int);

  // Description:
  // If true, then extract eigenvalues from tensor. False by default.
  vtkGetMacro(ExtractEigenvalues, int);
  vtkSetMacro(ExtractEigenvalues, int);
  vtkBooleanMacro(ExtractEigenvalues, int);

  // Description:
  // Turn on/off coloring of glyph with input scalar data or
  // eigenvalues. If false, or input scalar data not present, then the
  // scalars from the source object are passed through the filter.
  vtkSetMacro(ColorGlyphs,int);
  vtkGetMacro(ColorGlyphs,int);
  vtkBooleanMacro(ColorGlyphs,int);

//BTX
  enum
  {
      COLOR_BY_SCALARS,
      COLOR_BY_EIGENVALUES
  };
//ETX

  // Description:
  // Set the color mode to be used for the glyphs.  This can be set to
  // use the input scalars (default) or to use the eigenvalues at the
  // point.  If ThreeGlyphs is set and the eigenvalues are chosen for
  // coloring then each glyph is colored by the corresponding
  // eigenvalue and if not set the color corresponding to the largest
  // eigenvalue is chosen.  The recognized values are:
  // COLOR_BY_SCALARS = 0 (default)
  // COLOR_BY_EIGENVALUES = 1
  vtkSetClampMacro(ColorMode, int, COLOR_BY_SCALARS, COLOR_BY_EIGENVALUES);
  vtkGetMacro(ColorMode, int);
  void SetColorModeToScalars()
    {this->SetColorMode(COLOR_BY_SCALARS);};
  void SetColorModeToEigenvalues()
    {this->SetColorMode(COLOR_BY_EIGENVALUES);};

  // Description:
  // Set the scale factor of the superquadric. Default is 1.
  vtkSetMacro(ScaleFactor,double);

  // Description:
  // Turn on/off scalar clamping. If scalar clamping is on, the ivar
  // MaxScaleFactor is used to control the maximum scale factor. (This is
  // useful to prevent uncontrolled scaling near singularities.)
  vtkSetMacro(ClampScaling,int);
  vtkGetMacro(ClampScaling,int);
  vtkBooleanMacro(ClampScaling,int);

  // Description:
  // Set/Get the maximum allowable scale factor. This value is compared to the
  // combination of the scale factor times the eigenvalue. If less, the scale
  // factor is reset to the MaxScaleFactor. The boolean ClampScaling has to
  // be "on" for this to work.
  vtkSetMacro(MaxScaleFactor,double);
  vtkGetMacro(MaxScaleFactor,double);


protected:
  vtkSuperquadricTensorGlyphFilter();
  ~vtkSuperquadricTensorGlyphFilter();

  /* implementation of algorithm */
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  virtual int FillInputPortInformation(int port, vtkInformation* info) override;
  void SetActiveTensors(int, int, int, int, const char *);

  int ThetaResolution;
  int PhiResolution;
  double ThetaRoundness;
  double PhiRoundness;
  double ScaleFactor;
  double Gamma;
  int FixedThetaPhiRoundness;
  int ExtractEigenvalues;
  int ColorGlyphs; // Boolean controls coloring with input scalar data
  int ColorMode; // The coloring mode to use for the glyphs.
  int ClampScaling; // Boolean controls whether scaling is clamped.
  double MaxScaleFactor; // Maximum scale factor (ScaleFactor*eigenvalue)

private:
  vtkSuperquadricTensorGlyphFilter(const vtkSuperquadricTensorGlyphFilter&);  // Not implemented.
  void operator=(const vtkSuperquadricTensorGlyphFilter&);  // Not implemented.
};

#endif
