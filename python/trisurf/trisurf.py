#!/usr/bin/python3

import configobj
import xml.etree.ElementTree as ET
import base64
import zlib
import io


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
		self.config=configobj.ConfigObj(io.StringIO(string))
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

	def initFromSnapshot(self, snapshotfile):
		try:
			tree = ET.parse(snapshotfile)
		except:
			print("Error reading snapshot file")
			exit(1)

		root = tree.getroot()
		tapetxt=root.find('tape')
		version=root.find('trisurfversion')
		#print("Reading snapshot made from: "+version.text)
		self.tape=Tape()
		#print(tapetxt.text)
		self.tape.setTape(tapetxt.text)

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

