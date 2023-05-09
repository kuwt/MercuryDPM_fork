# trace generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
body_ = XMLUnstructuredGridReader(registrationName='body_*', FileName=['D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_0.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_1.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_2.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_3.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_4.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_5.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_6.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_7.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_8.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_9.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_10.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_11.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_12.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_13.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_14.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_15.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_16.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_17.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_18.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_19.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_20.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_21.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_22.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_23.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_24.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_25.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_26.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_27.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_28.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_29.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_30.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_31.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_32.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_33.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_34.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_35.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_36.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_37.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_38.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_39.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_40.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_41.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_42.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_43.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_44.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_45.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_46.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_47.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_48.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_49.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_50.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_51.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_52.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_53.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_54.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_55.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_56.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_57.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_58.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_59.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_60.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_61.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_62.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_63.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_64.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_65.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_66.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_67.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_68.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_69.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_70.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_71.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_72.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_73.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_74.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_75.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_76.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_77.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_78.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_79.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_80.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_81.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_82.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_83.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_84.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_85.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_86.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_87.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_88.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_89.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_90.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_91.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_92.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_93.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_94.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_95.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_96.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_97.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_98.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_99.vtu'])
body_.PointArrayStatus = ['Position', 'Velocity', 'Radius']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on body_
body_.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
body_Display = Show(body_, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
body_Display.Representation = 'Surface'
body_Display.ColorArrayName = [None, '']
body_Display.SelectTCoordArray = 'None'
body_Display.SelectNormalArray = 'None'
body_Display.SelectTangentArray = 'None'
body_Display.OSPRayScaleArray = 'Position'
body_Display.OSPRayScaleFunction = 'PiecewiseFunction'
body_Display.SelectOrientationVectors = 'None'
body_Display.ScaleFactor = 9.46886978149414
body_Display.SelectScaleArray = 'None'
body_Display.GlyphType = 'Arrow'
body_Display.GlyphTableIndexArray = 'None'
body_Display.GaussianRadius = 0.473443489074707
body_Display.SetScaleArray = ['POINTS', 'Position']
body_Display.ScaleTransferFunction = 'PiecewiseFunction'
body_Display.OpacityArray = ['POINTS', 'Position']
body_Display.OpacityTransferFunction = 'PiecewiseFunction'
body_Display.DataAxesGrid = 'GridAxesRepresentation'
body_Display.PolarAxes = 'PolarAxesRepresentation'
body_Display.ScalarOpacityUnitDistance = 142.57504821347788
body_Display.OpacityArrayName = ['POINTS', 'Position']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
body_Display.ScaleTransferFunction.Points = [-43.18429946899414, 0.0, 0.5, 0.0, 43.44380187988281, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
body_Display.OpacityTransferFunction.Points = [-43.18429946899414, 0.0, 0.5, 0.0, 43.44380187988281, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=body_,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 9.46886978149414
glyph1.GlyphTransform = 'Transform2'

# show data in view
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.SelectTCoordArray = 'None'
glyph1Display.SelectNormalArray = 'None'
glyph1Display.SelectTangentArray = 'None'
glyph1Display.OSPRayScaleArray = 'Position'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 9.572396850585937
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 0.4786198425292969
glyph1Display.SetScaleArray = ['POINTS', 'Position']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Position']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-43.18429946899414, 0.0, 0.5, 0.0, 43.07080078125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-43.18429946899414, 0.0, 0.5, 0.0, 43.07080078125, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1Display
glyph1Display.Specular = 1.0

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'
glyph1.ScaleArray = ['POINTS', 'Radius']
glyph1.ScaleFactor = 2.0
glyph1.GlyphMode = 'All Points'

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data bounds
renderView1.ResetCamera(-48.75773620605469, 48.07275390625, -31.664987564086914, 37.105674743652344, -49.08180236816406, 54.85300064086914)

# Properties modified on renderView1
renderView1.UseGradientBackground = 1

# change solid color
glyph1Display.AmbientColor = [1.0, 1.0, 0.4980392156862745]
glyph1Display.DiffuseColor = [1.0, 1.0, 0.4980392156862745]

animationScene1.Play()