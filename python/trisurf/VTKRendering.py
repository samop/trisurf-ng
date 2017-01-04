import os,sys
from . import trisurf
if sys.version_info<(3,0):
	from vtk import *



class MultiRender:
	def __init__(self,args,host):
		target_runs=getRargetRunIdxList(args)
		if target_runs==None:
			target_runs=list(range(1,len(host['runs'])+1))
		nruns=len(target_runs)
		#prepare rendering window
		self.renderer_window = vtkRenderWindow()
		self.renderer_window.AddRenderer(self.renderer)
		self.renderer_window.SetSize(1200,600)
		interactor = vtkRenderWindowInteractor()
		interactor.SetRenderWindow(self.renderer_window)
		interactor.Initialize()
		interactor.AddObserver("TimerEvent", self.RenderUpdate)
		timerIDR = interactor.CreateRepeatingTimer(1000)
		self.filename=[]
		self.renderer=[]
		i=0
		for run in target_runs:	
			#for each target run calculate renderer location
			r.vtkRenderer()
			r.SetBackground(0,0,0)
			p=1.0/float(nruns)
			x1=i*p
			x2=(i+1)*p
			r.SetViewport(x1,0.0,x2,1.0)
			self.renderer.Append(r)
			self.renderer_window.AddRenderer(r)
			i=i+1
		#call Renderer object with Run, renderer
		#start endless loop of interactor
		interactor.Start()

	def lastVTU(self,run):
		Dir=trisurf.Directory(maindir=run.maindir,simdir=run.subdir)
		filename=os.path.join("./",Dir.fullpath(),run.getLastVTU())
		return filename

	def lastActorForRun(self,run):
		filename=self.lastVTU(run)
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
		i=0
		for run in runs:
			if(self.lastVTU(run)!=self.filename[i]):
				#print("updejt")
				self.renderer.RemoveActor(self.actor)
				self.renderer.RemoveActor(self.textactor)
				self.actor=self.lastActor()
				self.textactor=self.textActor()
				self.renderer.AddActor(self.actor)
				self.renderer.AddActor(self.textactor)
				self.renderer_window.Render()
			#self.render.RemoveActor(self.actor)
			i=i+1	
		return


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
