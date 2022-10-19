# trace generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'STL Reader'
voxelsstl = STLReader(registrationName='voxels.stl', FileNames=['D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\voxelgrid\\voxels.stl'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
voxelsstlDisplay = Show(voxelsstl, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
voxelsstlDisplay.Representation = 'Surface'
voxelsstlDisplay.ColorArrayName = [None, '']
voxelsstlDisplay.SelectTCoordArray = 'None'
voxelsstlDisplay.SelectNormalArray = 'None'
voxelsstlDisplay.SelectTangentArray = 'None'
voxelsstlDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
voxelsstlDisplay.SelectOrientationVectors = 'None'
voxelsstlDisplay.ScaleFactor = 0.30000000000000004
voxelsstlDisplay.SelectScaleArray = 'None'
voxelsstlDisplay.GlyphType = 'Arrow'
voxelsstlDisplay.GlyphTableIndexArray = 'None'
voxelsstlDisplay.GaussianRadius = 0.015
voxelsstlDisplay.SetScaleArray = [None, '']
voxelsstlDisplay.ScaleTransferFunction = 'PiecewiseFunction'
voxelsstlDisplay.OpacityArray = [None, '']
voxelsstlDisplay.OpacityTransferFunction = 'PiecewiseFunction'
voxelsstlDisplay.DataAxesGrid = 'GridAxesRepresentation'
voxelsstlDisplay.PolarAxes = 'PolarAxesRepresentation'

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data bounds
renderView1.ResetCamera(-1.5, 1.5, -0.5099999904632568, 0.5099999904632568, -0.5099999904632568, 0.5099999904632568)

# Properties modified on renderView1
renderView1.UseGradientBackground = 1

# change solid color
voxelsstlDisplay.AmbientColor = [1.0, 0.0, 0.0]
voxelsstlDisplay.DiffuseColor = [1.0, 0.0, 0.0]

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(982, 745)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [-3.306197709934976, 4.310159039991645, 3.441757407325509]
renderView1.CameraViewUp = [0.29775249495107825, 0.7246527932304414, -0.6214674416361065]
renderView1.CameraParallelScale = 1.6643917749571595

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).