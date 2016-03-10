#!/usr/bin/python3

import configobj
import xml.etree.ElementTree as ET
import base64
import zlib
import io
import os
from itertools import islice
import mmap
import shlex
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

	def __str__(self):
		retval=""
		for key,val in self.config.iteritems():
			retval=retval+str(key)+" = "+str(val)+"\n"
		return retval



class Directory:
	def __init__(self, maindir=".", simdir=""):
		self.maindir=maindir
		self.simdir=simdir
		return

	def fullpath(self):
		return os.path.join(self.maindir,self.simdir)

	def exists(self):
		path=self.fullpath()
		if(os.path.exists(path)):
			return 1
		else:
			return 0

	def make(self):
		try:
			os.makedirs(self.fullpath())
		except:
			print("Cannot make directory "+self.fullpath()+"\n")
			exit(1)
		return

	def makeifnotexist(self):
		if(self.exists()==0):
			self.make()
		return

	def remove(self):
		if(self.exists()):
			try:
				os.rmdir(self.fullpath())
			except:
				print("Cannot remove directory "+self.fullpath()+ "\n")
				exit(1)
		return

	def goto(self):
		try:
			os.chdir(self.fullpath())
		except:
			print("Cannot go to directory "+self.fullpath()+"\n")
		return


class Statistics:
	def __init__(self,path,filename="statistics.csv"):
		self.path=path
		self.filename=filename
		self.fullname=os.path.join(path,filename)
		self.read()
		return

	def __str__(self):
		return(str(self.fullname))

	def exists(self):
		if(os.path.isfile(self.fullname)):
			return True
		else:
			return False

	def mapcount(self):
		f = open(self.fullname, "r+")
		buf = mmap.mmap(f.fileno(), 0)
		lines = 0
		readline = buf.readline
		while readline():
			lines += 1
		return lines

	def read(self):
		if(self.exists()):
			nlines=self.mapcount()
			try:
				with open(self.fullname, "r+") as fin:
					i=0;
					for line in fin:
						if(i==1):
							#print (line)
							fields=shlex.split(line)
							epoch1=fields[0]
							n1=fields[1]
						if(i==nlines-1):
							fields=shlex.split(line)
							epoch2=fields[0]
							n2=fields[1]
						i=i+1
			except:
				print("Cannot read statistics file in "+self.fullname+"\n")
				return(False)
		else:
			print("File "+self.fullname+" does not exists.\n")
			return(False)

		self.dT=(int(epoch2)-int(epoch1))/(int(n2)-int(n1))
		return(True)


class Runner:
	'''
	Class Runner consists of a single running or terminated instance of the trisurf
	'''
	def __init__(self, subdir='run0', tape='', snapshot=''):
		self.subdir=subdir
		if(tape!=''):
			self.initFromTape(tape)
		if(snapshot!=''):
			self.initFromSnapshot(snapshot)
		return


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
		self.tape=Tape()
		self.tape.setTape(tapetxt.text)

	def getStatus(self):
		return 0

	def start(self):
		if(self.getStatus()==0):
			self.Dir=Directory(maindir=self.maindir,simdir=self.subdir)
			self.Dir.makeifnotexist()
			self.Dir.goto()
			print("Starting trisurf-ng executable at "+self.Dir.fullpath()+"\n")
		else:
			print("Process already running. Not starting\n")
		return

	def stop(self):
		pass

	def setMaindir(self,prefix,variables):
		maindir="./"
		for p,v in zip(prefix,variables):
			if(v=="xk0"):
				tv=str(round(float(self.tape.config[v])))
			else:
				tv=self.tape.config[v]
			maindir=maindir+p+tv
		self.maindir=maindir
		return

	def setSubdir(self, subdir="run0"):
		self.subdir=subdir
		return

	def getStatistics(self, statfile="statistics.csv"):
		self.statistics=Statistics("", statfile) # we are already in the running directory, so local path is needed!
		return

	def __str__(self):
		if(self.getStatus()==0):
			str=" not running."
		else:
			str=" running."
		return(self.Dir.fullpath()+str)

