#script to visualise the output of FE runs in paraview.
# usage: pvbatch --use-offscreen-rendering  makeMovie.py
from paraview.simple import *
import os
import glob

#mkdir -p vtk && mv *Processor_[6-9]*.vtu vtk  && mv *Processor_[3-5]*.vtu vtk  && mv *Processor_[1-2]*.vtu vtk  && mv *Processor_[0-2]*.vtu vtk

# get active view
renderView = GetActiveViewOrCreate('RenderView')
# make background white
renderView.Background = [1.0, 1.0, 1.0]
# Turn camera into sideway view
#renderView.ResetCamera(-0.00249999994412, 0.496250003576, -0.035000000149, 0.0960000008345, -0.102499999106, 0.354999989271)
#camera=GetActiveCamera()
#camera.Elevation(-90)
renderView.InteractionMode = '2D'
renderView.CameraPosition = [0.24750000180210918, -1.8466974306691022, 0.1262499950826168]
renderView.CameraFocalPoint = [0.24750000180210918, 0.029577700421214104, 0.1262499950826168]
renderView.CameraViewUp = [0.0, 0.0, 1.0]

ids = ['1','2','2b','3','4','5','6','8','11','14','17','19','19b','21','21b','22','22b','30','82']
#ids = ['14']

for id in ids: 
	print 'id: ',id

	movieDir = 'v9movie'
	vtkDir = 'v9run'+id
	
	#make jpg directory
	try:
		os.mkdir(movieDir)
	except OSError:
		print('')


	#Find the maximum timestep
	Data = glob.glob(vtkDir+'/GCG'+id+'Particle_*.vtu')
	maxTime = 0
	for fileName in Data:
		tokens1 = fileName.split('.')
		tokens2 = tokens1[-2].split('_')
		if int(tokens2[-1]) > maxTime:
			maxTime = int(tokens2[-1])
	#maxTime = 200
	print '#timesteps ', maxTime

	# snapshot every n-th timestep
	step=10
	for i in range(0,maxTime+1,step):
		particleFile = vtkDir+'/GCG'+id+'Particle_'+str(i)+'.vtu'
		wallFile = vtkDir+'/GCG'+id+'Wall_'+str(i)+'.vtu'
		imageFile = movieDir+'/GCG'+id+'_'+str(i)+'.png'
		
		if len(glob.glob(imageFile)):
			#print 'snapshot:', i, 'already exists'
			continue
		
		print 'snapshot:', i
		#particles
		particles = XMLUnstructuredGridReader(FileName=particleFile)
		glyph = Glyph(particles)
		glyph.GlyphType = 'Sphere'
		glyph.Scalars = 'Radius'
		glyph.Vectors = 'None'
		glyph.ScaleMode = 'scalar'
		glyph.ScaleFactor = 2
		glyph.GlyphMode = 'All Points'
		Show(glyph)	
		#Walls
		DataW = glob.glob(wallFile)
		walls = XMLUnstructuredGridReader(FileName=DataW)
		Show(walls)
		# show data in view
		#text = Text(Text=id)
		#textDisplay = GetDisplayProperties(text)
		#textDisplay.Color = [0.0, 0.0, 0.0]
		#textDisplay.WindowLocation = 'LowerCenter'
		#Show(text)

		# Render
		Render()
		#Interact()
		SaveScreenshot(imageFile, magnification=2, view=renderView)

		#delete pv data (not sure if this is done right)
		Delete(particles)
		del particles
		Delete(glyph)
		del glyph
		Delete(walls)
		del walls

	#make movie
	imageFiles = movieDir+'/GCG'+id+'_%d00.png'
	movieFile = movieDir+'/GCG'+id+'.avi'
	if len(glob.glob(imageFile)):
		print 'video:', movieFile, 'already exists'
		continue
	os.system('ffmpeg -y -i '+imageFiles+' -c:v libxvid -qscale:v 9 '+movieFile) 
#finally, transfer it to home: 
os.system('rsync -avz v9movie/ ut148149.roaming.utwente.nl:MercuryLab/Presentations/Movies/v9movie')

