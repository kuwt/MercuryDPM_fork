/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTensorGlyphSameEigensystem.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTensorGlyphSameEigensystem - scale and orient a single glyph according to tensor eigenvalues and eigenvectors
// .SECTION Description
// vtkTensorGlyphSameEigensystem is a filter that copies a geometric representation
// (specified as polygonal data) to all its input points and uses the solution
// of the same eigenvalue problem to orient and scale them. The differences
// compared to its superclass vtkTensorGlyph are threefold: First, the fact
// that only a single eigenvalue problem is involved, not a eigenvalue problem
// in every point. Second, the possibility to pass the single solution of the
// eigenvalue equation and skip solving it. Third, the fact that it names the
// data array it creates such that ParaView will be able to use them directly
// for coloring and further processing.

// .SECTION Thanks
// Class has been developed and contributed by Sven Buijssen, TU Dortmund,
// Germany.

// .SECTION See Also
// vtkTensorGlyph vtkGlyph3D vtkPointLoad vtkHyperStreamline

#ifndef __vtkTensorGlyphSameEigensystem_h
#define __vtkTensorGlyphSameEigensystem_h

#include "vtkTensorGlyph.h"

// Windows specific stuff------------------------------------------
#if defined(_WIN32) || defined(WIN32) || defined(__CYGWIN__)
#  if defined(VTK_BUILD_SHARED_LIBS)
#    define VTK_SUPERQUADRICTENSORGLYPH_EXPORT __declspec( dllexport )
#  else
#    define VTK_SUPERQUADRICTENSORGLYPH_EXPORT
#  endif
#else
#  define VTK_SUPERQUADRICTENSORGLYPH_EXPORT
#endif


class VTK_SUPERQUADRICTENSORGLYPH_EXPORT vtkTensorGlyphSameEigensystem : public vtkTensorGlyph
{
public:
  vtkTypeMacro(vtkTensorGlyphSameEigensystem,vtkTensorGlyph);

  // Description
  // Construct object with scaling on and scale factor 1.0. Eigenvalues are
  // extracted, glyphs are colored with input scalar data, and logarithmic
  // scaling is turned off.
  static vtkTensorGlyphSameEigensystem *New();

  void SetEigenvectorsEigenValues(double **ev, double *ew);

protected:
  vtkTensorGlyphSameEigensystem();
  ~vtkTensorGlyphSameEigensystem();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
//  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

//BTX
  int HaveEigenvalues; // Boolean controls whether eigenfunction extraction is required
  double Eigenvalues[3];
  double EigenvectorsX[3], EigenvectorsY[3], EigenvectorsZ[3];
//ETX


private:
  vtkTensorGlyphSameEigensystem(const vtkTensorGlyphSameEigensystem&);  // Not implemented.
  void operator=(const vtkTensorGlyphSameEigensystem&);  // Not implemented.
};

#endif
