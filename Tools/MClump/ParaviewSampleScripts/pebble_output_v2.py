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

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'
glyph1.ScaleArray = ['POINTS', 'Radius']
glyph1.ScaleFactor = 2.0

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data bounds
renderView1.ResetCamera(-48.75773620605469, 46.30080795288086, -31.658517837524414, 34.493350982666016, -49.08180236816406, 54.80889892578125)

# Properties modified on glyph1
glyph1.GlyphMode = 'All Points'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1Display
glyph1Display.Specular = 1.0

# Properties modified on renderView1
renderView1.UseGradientBackground = 1

# set active source
SetActiveSource(None)

# create a new 'XML Unstructured Grid Reader'
body_pd_ = XMLUnstructuredGridReader(registrationName='body_pd_*', FileName=['D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_0.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_1.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_2.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_3.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_4.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_5.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_6.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_7.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_8.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_9.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_10.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_11.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_12.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_13.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_14.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_15.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_16.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_17.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_18.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_19.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_20.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_21.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_22.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_23.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_24.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_25.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_26.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_27.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_28.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_29.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_30.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_31.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_32.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_33.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_34.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_35.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_36.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_37.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_38.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_39.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_40.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_41.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_42.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_43.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_44.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_45.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_46.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_47.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_48.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_49.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_50.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_51.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_52.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_53.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_54.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_55.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_56.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_57.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_58.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_59.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_60.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_61.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_62.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_63.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_64.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_65.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_66.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_67.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_68.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_69.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_70.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_71.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_72.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_73.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_74.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_75.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_76.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_77.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_78.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_79.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_80.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_81.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_82.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_83.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_84.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_85.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_86.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_87.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_88.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_89.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_90.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_91.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_92.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_93.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_94.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_95.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_96.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_97.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_98.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_99.vtu'])
body_pd_.PointArrayStatus = ['Position', 'v1', 'v2', 'v3', 'Radius']

# Properties modified on body_pd_
body_pd_.TimeArray = 'None'

# show data in view
body_pd_Display = Show(body_pd_, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
body_pd_Display.Representation = 'Surface'
body_pd_Display.ColorArrayName = [None, '']
body_pd_Display.SelectTCoordArray = 'None'
body_pd_Display.SelectNormalArray = 'None'
body_pd_Display.SelectTangentArray = 'None'
body_pd_Display.OSPRayScaleArray = 'Position'
body_pd_Display.OSPRayScaleFunction = 'PiecewiseFunction'
body_pd_Display.SelectOrientationVectors = 'None'
body_pd_Display.ScaleFactor = 0.1
body_pd_Display.SelectScaleArray = 'None'
body_pd_Display.GlyphType = 'Arrow'
body_pd_Display.GlyphTableIndexArray = 'None'
body_pd_Display.GaussianRadius = 0.005
body_pd_Display.SetScaleArray = ['POINTS', 'Position']
body_pd_Display.ScaleTransferFunction = 'PiecewiseFunction'
body_pd_Display.OpacityArray = ['POINTS', 'Position']
body_pd_Display.OpacityTransferFunction = 'PiecewiseFunction'
body_pd_Display.DataAxesGrid = 'GridAxesRepresentation'
body_pd_Display.PolarAxes = 'PolarAxesRepresentation'
body_pd_Display.ScalarOpacityUnitDistance = 0.0
body_pd_Display.OpacityArrayName = ['POINTS', 'Position']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
body_pd_Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
body_pd_Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Glyph'
glyph2 = Glyph(registrationName='Glyph2', Input=body_pd_,
    GlyphType='Arrow')
glyph2.OrientationArray = ['POINTS', 'No orientation array']
glyph2.ScaleArray = ['POINTS', 'No scale array']
glyph2.ScaleFactor = 0.1
glyph2.GlyphTransform = 'Transform2'

# Properties modified on glyph2
glyph2.OrientationArray = ['POINTS', 'v1']
glyph2.ScaleFactor = 100.0

# show data in view
glyph2Display = Show(glyph2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph2Display.Representation = 'Surface'
glyph2Display.ColorArrayName = [None, '']
glyph2Display.SelectTCoordArray = 'None'
glyph2Display.SelectNormalArray = 'None'
glyph2Display.SelectTangentArray = 'None'
glyph2Display.OSPRayScaleArray = 'Position'
glyph2Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph2Display.SelectOrientationVectors = 'None'
glyph2Display.ScaleFactor = 10.0
glyph2Display.SelectScaleArray = 'None'
glyph2Display.GlyphType = 'Arrow'
glyph2Display.GlyphTableIndexArray = 'None'
glyph2Display.GaussianRadius = 0.5
glyph2Display.SetScaleArray = ['POINTS', 'Position']
glyph2Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph2Display.OpacityArray = ['POINTS', 'Position']
glyph2Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph2Display.DataAxesGrid = 'GridAxesRepresentation'
glyph2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change solid color
glyph2Display.AmbientColor = [1.0, 0.0, 0.0]
glyph2Display.DiffuseColor = [1.0, 0.0, 0.0]

# Properties modified on glyph2Display
glyph2Display.Specular = 1.0

# set active source
SetActiveSource(None)

# create a new 'XML Unstructured Grid Reader'
body_pd__1 = XMLUnstructuredGridReader(registrationName='body_pd_*', FileName=['D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_0.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_1.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_2.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_3.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_4.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_5.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_6.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_7.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_8.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_9.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_10.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_11.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_12.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_13.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_14.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_15.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_16.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_17.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_18.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_19.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_20.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_21.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_22.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_23.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_24.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_25.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_26.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_27.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_28.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_29.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_30.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_31.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_32.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_33.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_34.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_35.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_36.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_37.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_38.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_39.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_40.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_41.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_42.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_43.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_44.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_45.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_46.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_47.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_48.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_49.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_50.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_51.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_52.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_53.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_54.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_55.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_56.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_57.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_58.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_59.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_60.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_61.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_62.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_63.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_64.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_65.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_66.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_67.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_68.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_69.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_70.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_71.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_72.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_73.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_74.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_75.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_76.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_77.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_78.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_79.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_80.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_81.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_82.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_83.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_84.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_85.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_86.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_87.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_88.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_89.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_90.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_91.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_92.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_93.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_94.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_95.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_96.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_97.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_98.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_99.vtu'])
body_pd__1.PointArrayStatus = ['Position', 'v1', 'v2', 'v3', 'Radius']

# Properties modified on body_pd__1
body_pd__1.TimeArray = 'None'

# show data in view
body_pd__1Display = Show(body_pd__1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
body_pd__1Display.Representation = 'Surface'
body_pd__1Display.ColorArrayName = [None, '']
body_pd__1Display.SelectTCoordArray = 'None'
body_pd__1Display.SelectNormalArray = 'None'
body_pd__1Display.SelectTangentArray = 'None'
body_pd__1Display.OSPRayScaleArray = 'Position'
body_pd__1Display.OSPRayScaleFunction = 'PiecewiseFunction'
body_pd__1Display.SelectOrientationVectors = 'None'
body_pd__1Display.ScaleFactor = 0.1
body_pd__1Display.SelectScaleArray = 'None'
body_pd__1Display.GlyphType = 'Arrow'
body_pd__1Display.GlyphTableIndexArray = 'None'
body_pd__1Display.GaussianRadius = 0.005
body_pd__1Display.SetScaleArray = ['POINTS', 'Position']
body_pd__1Display.ScaleTransferFunction = 'PiecewiseFunction'
body_pd__1Display.OpacityArray = ['POINTS', 'Position']
body_pd__1Display.OpacityTransferFunction = 'PiecewiseFunction'
body_pd__1Display.DataAxesGrid = 'GridAxesRepresentation'
body_pd__1Display.PolarAxes = 'PolarAxesRepresentation'
body_pd__1Display.ScalarOpacityUnitDistance = 0.0
body_pd__1Display.OpacityArrayName = ['POINTS', 'Position']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
body_pd__1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
body_pd__1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Glyph'
glyph3 = Glyph(registrationName='Glyph3', Input=body_pd__1,
    GlyphType='Arrow')
glyph3.OrientationArray = ['POINTS', 'No orientation array']
glyph3.ScaleArray = ['POINTS', 'No scale array']
glyph3.ScaleFactor = 0.1
glyph3.GlyphTransform = 'Transform2'

# Properties modified on glyph3
glyph3.OrientationArray = ['POINTS', 'v2']
glyph3.ScaleFactor = 100.0
glyph3.GlyphMode = 'All Points'

# show data in view
glyph3Display = Show(glyph3, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph3Display.Representation = 'Surface'
glyph3Display.ColorArrayName = [None, '']
glyph3Display.SelectTCoordArray = 'None'
glyph3Display.SelectNormalArray = 'None'
glyph3Display.SelectTangentArray = 'None'
glyph3Display.OSPRayScaleArray = 'Position'
glyph3Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph3Display.SelectOrientationVectors = 'None'
glyph3Display.ScaleFactor = 10.0
glyph3Display.SelectScaleArray = 'None'
glyph3Display.GlyphType = 'Arrow'
glyph3Display.GlyphTableIndexArray = 'None'
glyph3Display.GaussianRadius = 0.5
glyph3Display.SetScaleArray = ['POINTS', 'Position']
glyph3Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph3Display.OpacityArray = ['POINTS', 'Position']
glyph3Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph3Display.DataAxesGrid = 'GridAxesRepresentation'
glyph3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change solid color
glyph3Display.AmbientColor = [1.0, 1.0, 0.0]
glyph3Display.DiffuseColor = [1.0, 1.0, 0.0]

# set active source
SetActiveSource(None)

# create a new 'XML Unstructured Grid Reader'
body_pd__2 = XMLUnstructuredGridReader(registrationName='body_pd_*', FileName=['D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_0.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_1.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_2.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_3.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_4.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_5.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_6.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_7.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_8.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_9.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_10.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_11.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_12.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_13.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_14.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_15.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_16.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_17.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_18.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_19.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_20.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_21.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_22.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_23.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_24.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_25.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_26.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_27.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_28.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_29.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_30.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_31.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_32.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_33.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_34.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_35.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_36.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_37.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_38.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_39.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_40.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_41.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_42.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_43.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_44.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_45.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_46.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_47.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_48.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_49.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_50.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_51.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_52.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_53.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_54.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_55.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_56.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_57.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_58.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_59.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_60.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_61.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_62.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_63.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_64.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_65.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_66.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_67.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_68.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_69.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_70.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_71.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_72.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_73.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_74.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_75.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_76.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_77.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_78.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_79.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_80.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_81.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_82.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_83.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_84.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_85.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_86.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_87.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_88.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_89.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_90.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_91.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_92.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_93.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_94.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_95.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_96.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_97.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_98.vtu', 'D:\\MDPM_Multiparticle\\MercurySource\\Tools\\MClump\\output\\paraview\\body_pd_99.vtu'])
body_pd__2.PointArrayStatus = ['Position', 'v1', 'v2', 'v3', 'Radius']

# Properties modified on body_pd__2
body_pd__2.TimeArray = 'None'

# show data in view
body_pd__2Display = Show(body_pd__2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
body_pd__2Display.Representation = 'Surface'
body_pd__2Display.ColorArrayName = [None, '']
body_pd__2Display.SelectTCoordArray = 'None'
body_pd__2Display.SelectNormalArray = 'None'
body_pd__2Display.SelectTangentArray = 'None'
body_pd__2Display.OSPRayScaleArray = 'Position'
body_pd__2Display.OSPRayScaleFunction = 'PiecewiseFunction'
body_pd__2Display.SelectOrientationVectors = 'None'
body_pd__2Display.ScaleFactor = 0.1
body_pd__2Display.SelectScaleArray = 'None'
body_pd__2Display.GlyphType = 'Arrow'
body_pd__2Display.GlyphTableIndexArray = 'None'
body_pd__2Display.GaussianRadius = 0.005
body_pd__2Display.SetScaleArray = ['POINTS', 'Position']
body_pd__2Display.ScaleTransferFunction = 'PiecewiseFunction'
body_pd__2Display.OpacityArray = ['POINTS', 'Position']
body_pd__2Display.OpacityTransferFunction = 'PiecewiseFunction'
body_pd__2Display.DataAxesGrid = 'GridAxesRepresentation'
body_pd__2Display.PolarAxes = 'PolarAxesRepresentation'
body_pd__2Display.ScalarOpacityUnitDistance = 0.0
body_pd__2Display.OpacityArrayName = ['POINTS', 'Position']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
body_pd__2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
body_pd__2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Glyph'
glyph4 = Glyph(registrationName='Glyph4', Input=body_pd__2,
    GlyphType='Arrow')
glyph4.OrientationArray = ['POINTS', 'No orientation array']
glyph4.ScaleArray = ['POINTS', 'No scale array']
glyph4.ScaleFactor = 0.1
glyph4.GlyphTransform = 'Transform2'

# Properties modified on glyph4
glyph4.OrientationArray = ['POINTS', 'v3']
glyph4.ScaleFactor = 100.0
glyph4.GlyphMode = 'All Points'

# show data in view
glyph4Display = Show(glyph4, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph4Display.Representation = 'Surface'
glyph4Display.ColorArrayName = [None, '']
glyph4Display.SelectTCoordArray = 'None'
glyph4Display.SelectNormalArray = 'None'
glyph4Display.SelectTangentArray = 'None'
glyph4Display.OSPRayScaleArray = 'Position'
glyph4Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph4Display.SelectOrientationVectors = 'None'
glyph4Display.ScaleFactor = 10.0
glyph4Display.SelectScaleArray = 'None'
glyph4Display.GlyphType = 'Arrow'
glyph4Display.GlyphTableIndexArray = 'None'
glyph4Display.GaussianRadius = 0.5
glyph4Display.SetScaleArray = ['POINTS', 'Position']
glyph4Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph4Display.OpacityArray = ['POINTS', 'Position']
glyph4Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph4Display.DataAxesGrid = 'GridAxesRepresentation'
glyph4Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph4Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph4Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# change solid color
glyph4Display.AmbientColor = [0.3333333333333333, 1.0, 0.0]
glyph4Display.DiffuseColor = [0.3333333333333333, 1.0, 0.0]

# set active source
SetActiveSource(glyph1)

# change solid color
glyph1Display.AmbientColor = [1.0, 0.6666666666666666, 0.0]
glyph1Display.DiffuseColor = [1.0, 0.6666666666666666, 0.0]

# reset view to fit data bounds
renderView1.ResetCamera(-48.75773620605469, 48.07275390625, -31.664987564086914, 37.105674743652344, -49.08180236816406, 54.85300064086914)