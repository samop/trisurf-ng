#!/usr/bin/python3

import configobj

'''
This is a trisurf instance manager written in python


Invoke with:
tsmgr [-t tape | -r snapshot.vtu] [-s subdirectory]

If tape is specified, the trisurf wilt start from tape with initial distribution, if snapshot is specified the trisurf will be restored from given snapshot file and simulation will continue.

'''



class Tape:
	'''Has all the info on the tape'''

	def __init__(self):
		return

	def readTape(self, tape='tape'):
		try:
			self.config=configobj.ConfigObj(tape)
		except:
			print("Error reading or parsing tape file!\n")
			exit(1)
	

	def setTape(self, string):
		self.tape=string
		return

	def getValue(self,key):
		return self.config[key]



class Runner:
	'''
	Class Runner consists of a single running or terminated instance of the trisurf
	'''
	def initFromTape(self, tape):
		self.tape=Tape()
		self.tape.readTape(tape)
		pass

	def initFromSnapshot(self, tape='snapshot.vtu'):
		pass

	def __init__(self, subdir='run0', tape='', snapshot=''):
		self.subdir=subdir
		if(tape!=''):
			self.initFromTape(tape)
		if(snapshot!=''):
			self.initFromSnapshot(snapshot)

		return

	def getStatus(self):
		pass

	def start(self):
		pass

	def stop(self):
		pass

	def __str__(self):
		return("Running instance")

