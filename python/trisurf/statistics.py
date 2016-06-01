from trisurf import trisurf
import os

def combine(Runs):
	"""Runs are those runs of which statistics are to be combined"""
	data=[]
	for run in Runs:
		dir=trisurf.Directory(maindir=run.maindir,simdir=run.subdir)
		statfile=os.path.join(dir.fullpath(),"statistics.csv")
		with open (statfile,"r") as myfile:
			#lines = [line.rstrip('\n') for line in myfile]
			data=data+myfile.readlines()[1:]
	
	print("Header line placer... Not yet implemented")
	for line in data:
		print(line.replace('\n',''))

	

