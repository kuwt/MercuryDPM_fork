#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
drumSizeWall_0vtu = XMLUnstructuredGridReader(FileName=['/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeWall_0.vtu'])

# set active source
SetActiveSource(drumSizeWall_0vtu)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1018, 597]

# show data in view
drumSizeWall_0vtuDisplay = Show(drumSizeWall_0vtu, renderView1)
# trace defaults for the display properties.
drumSizeWall_0vtuDisplay.ColorArrayName = [None, '']
drumSizeWall_0vtuDisplay.ScalarOpacityUnitDistance = 0.6403124332846835

# reset view to fit data
renderView1.ResetCamera()

# create a new 'XML Unstructured Grid Reader'
drumSizeParticle_ = XMLUnstructuredGridReader(FileName=['/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_0.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_1.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_2.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_3.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_4.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_5.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_6.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_7.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_8.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_9.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_10.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_11.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_12.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_13.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_14.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_15.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_16.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_17.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_18.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_19.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_20.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_21.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_22.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_23.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_24.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_25.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_26.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_27.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_28.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_29.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_30.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_31.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_32.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_33.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_34.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_35.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_36.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_37.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_38.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_39.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_40.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_41.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_42.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_43.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_44.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_45.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_46.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_47.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_48.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_49.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_50.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_51.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_52.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_53.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_54.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_55.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_56.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_57.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_58.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_59.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_60.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_61.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_62.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_63.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_64.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_65.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_66.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_67.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_68.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_69.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_70.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_71.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_72.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_73.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_74.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_75.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_76.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_77.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_78.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_79.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_80.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_81.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_82.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_83.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_84.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_85.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_86.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_87.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_88.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_89.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_90.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_91.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_92.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_93.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_94.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_95.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_96.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_97.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_98.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_99.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_100.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_101.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_102.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_103.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_104.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_105.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_106.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_107.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_108.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_109.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_110.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_111.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_112.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_113.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_114.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_115.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_116.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_117.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_118.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_119.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_120.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_121.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_122.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_123.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_124.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_125.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_126.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_127.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_128.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_129.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_130.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_131.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_132.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_133.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_134.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_135.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_136.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_137.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_138.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_139.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_140.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_141.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_142.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_143.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_144.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_145.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_146.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_147.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_148.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_149.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_150.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_151.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_152.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_153.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_154.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_155.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_156.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_157.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_158.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_159.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_160.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_161.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_162.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_163.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_164.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_165.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_166.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_167.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_168.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_169.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_170.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_171.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_172.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_173.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_174.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_175.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_176.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_177.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_178.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_179.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_180.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_181.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_182.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_183.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_184.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_185.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_186.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_187.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_188.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_189.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_190.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_191.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_192.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_193.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_194.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_195.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_196.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_197.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_198.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_199.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_200.vtu', '/Users/antRthorn/Documents/UniWork/svn/code/MercuryDPM/TrunkBuild/Drivers/USER/CuriousU/DrumSizeParticle_201.vtu'])
drumSizeParticle_.PointArrayStatus = ['Velocity', 'Diameter', 'SpeciesType']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# show data in view
drumSizeParticle_Display = Show(drumSizeParticle_, renderView1)
# trace defaults for the display properties.
drumSizeParticle_Display.ColorArrayName = [None, '']
drumSizeParticle_Display.ScalarOpacityUnitDistance = 0.5042527701066964

# show data in view
drumSizeWall_0vtuDisplay = Show(drumSizeWall_0vtu, renderView1)

# create a new 'Glyph'
glyph1 = Glyph(Input=drumSizeParticle_,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'Diameter']
glyph1.Vectors = ['POINTS', 'Velocity']
glyph1.ScaleFactor = 0.036463400721549986
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = 2.0

# get color transfer function/color map for 'Diameter'
diameterLUT = GetColorTransferFunction('Diameter')

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.ColorArrayName = ['POINTS', 'Diameter']
glyph1Display.LookupTable = diameterLUT

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'Diameter'
diameterPWF = GetOpacityTransferFunction('Diameter')

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

animationScene1.Play()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.0, -1.2869886331524016, 0.0]
renderView1.CameraFocalPoint = [0.0, -0.05000000074505806, 0.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 0.32015621664234173

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
