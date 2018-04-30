# Usage:
#
#   Pipeline:
#     vtkPointData->SetTensors(data)
#     where data is a vtkDataArray of, e.g., type double, with 9 components
#     (the three rows of the 3x3 tensor)
#
#   Note:
#     derived from TensorGlyphFilter by Bryn Lloyd, blloyd@vision.ee.ethz.ch
#     (see http://markmail.org/message/cdrnxa2giayqfkvr) and the server
#     manager XML posted on the ParaView mailing list by Jean M. Favre,
#     jfavre@cscs.ch (see http://markmail.org/message/o5lt4er3dlxjpqqy)
#     and inspired by the paper
#       Gordon Kindlmann. Superquadric Tensor Glyphs.
#       In Proceedings IEEE TVCG/EG Symposium on Visualization 2004,
#         pages 147-154, May 2004, http://www.cs.utah.edu/~gk/papers/vissym04/
#
#
# Contribution :
#   Author:
#     Sven H.M. Buijssen, sven.buijssen@tu-dortmund.de
#
#   Affiliation:
#     Institute of Applied Mathematics, Faculty of Mathematics
#     TU Dortmund University
#     Vogelpothsweg 87, 44 227 Dortmund, Germany.
#
#   Tweaks/performance optimization:
#     Joseph Hennessey, via upload to ParaView Wiki

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIRS}
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

IF(PARAVIEW_BUILD_QT_GUI)
  # Extend the auto-generated panel
  QT5_WRAP_CPP(MOC_SRCS pqSuperquadricTensorGlyphPanel.h)
  #ADD_PARAVIEW_OBJECT_PANEL(IFACES IFACE_SRCS
  #  CLASS_NAME pqSuperquadricTensorGlyphPanel
  #  XML_NAME SuperquadricTensorGlyphFilter # name of SourceProxy in server manager XML
  #  XML_GROUP filters
  #  )

  ADD_PARAVIEW_PLUGIN(SuperquadricTensorGlyphFilter "3.0.1"
    REQUIRED_ON_SERVER
    REQUIRED_ON_CLIENT
    SERVER_MANAGER_XML vtkSuperquadricTensorGlyphFilter_SM.xml
    SERVER_MANAGER_SOURCES
        vtkSuperquadricTensorGlyphFilter.cxx vtkTensorGlyphSameEigensystem.cxx
    #GUI_INTERFACES ${IFACES}
    #GUI_SOURCES ${MOC_SRCS} ${UI_SRCS} ${IFACE_SRCS}
    #    pqSuperquadricTensorGlyphPanel.cxx
    )
ELSE(PARAVIEW_BUILD_QT_GUI)
  ADD_PARAVIEW_PLUGIN(SuperquadricTensorGlyphFilter "3.0.1"
    SERVER_MANAGER_SOURCES
        vtkSuperquadricTensorGlyphFilter.cxx vtkTensorGlyphSameEigensystem.cxx
    SERVER_MANAGER_XML vtkSuperquadricTensorGlyphFilter_SM.xml
    )
ENDIF(PARAVIEW_BUILD_QT_GUI)