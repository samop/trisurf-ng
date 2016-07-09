import os,sys
from . import trisurf
if sys.version_info<(3,0):
	from vtk import *


class Renderer:
	def __init__(self,args,host):
		self.host=host
		self.args=args
		self.renderer = vtkRenderer()
		self.actor=self.lastActor()
		self.textactor=self.textActor()
		self.renderer.AddActor(self.actor)
		self.renderer.AddActor(self.textactor)
		self.renderer.SetBackground(0, 0, 0) # Set background to white
		# Create the RendererWindow
		self.renderer_window = vtkRenderWindow()
		self.renderer_window.AddRenderer(self.renderer)
		self.renderer_window.SetSize(1200,600)
	
		self.renderer.SetViewport(0.0,0.0,0.5,1.0)
		rend=vtk.vtkRenderer()
		rend.AddActor(self.actor)
		rend.SetViewport(0.5,0.0,1.0,1.0)
		self.renderer_window.AddRenderer(rend)	
# Set up a check for aborting rendering.
		# Create the RendererWindowInteractor and display the vtk_file
		interactor = vtkRenderWindowInteractor()
		interactor.SetRenderWindow(self.renderer_window)
		interactor.Initialize()
 		interactor.AddObserver("TimerEvent", self.RenderUpdate)
		timerIDR = interactor.CreateRepeatingTimer(1000)
		interactor.Start()

		return

	def lastVTU(self):
		Dir=trisurf.Directory(maindir=self.host['runs'][0].maindir,simdir=self.host['runs'][0].subdir)
		filename=os.path.join("./",Dir.fullpath(),self.host['runs'][0].getLastVTU())
		return filename

	def textActor(self):
		textactor=vtkTextActor()
		textactor.SetInput(self.filename)
		tp=textactor.GetTextProperty()
		tp.SetColor(1,1,1)
		tp.SetFontSize(11)
		textactor.SetDisplayPosition(20,30)
		return textactor

	def lastActor(self):
		self.filename=self.lastVTU()
		reader=vtkXMLUnstructuredGridReader()
		reader.SetFileName(self.filename)
		reader.Update() # Needed because of GetScalarRange
		output = reader.GetOutput()
		scalar_range = output.GetScalarRange()
		mapper = vtkDataSetMapper()
		mapper.SetInput(output)
		mapper.SetScalarRange(scalar_range)

		# Create the Actor
		actor = vtkActor()
		actor.SetMapper(mapper)
		return actor


	def RenderUpdate(self, obj, event):
		if(self.lastVTU()!=self.filename):
			#print("updejt")
			self.renderer.RemoveActor(self.actor)
			self.renderer.RemoveActor(self.textactor)
			self.actor=self.lastActor()
			self.textactor=self.textActor()
			self.renderer.AddActor(self.actor)
			self.renderer.AddActor(self.textactor)
			self.renderer_window.Render()
		#self.render.RemoveActor(self.actor)
		
		return
