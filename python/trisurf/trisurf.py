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
import psutil
import time
import datetime
from subprocess import call, Popen

# Process status
TS_NOLOCK=0 # lock file does not exist
TS_NONEXISTANT=0 # process is not in the list of processes
TS_STOPPED=1 # the process is listed, but is in stopped state
TS_RUNNING=2 # process is running

class FileContent:
	'''
	Class is helpful for reading and writting the specific files.
	'''
	def __init__(self,filename):
		''' The instance is done by calling constructor FileContent(filename)

		The object then reads filename if it exists, otherwise the data is empty string.
		User may force to reread file by calling the readline() method of the class.

		Filename is stored in local variable for future file operations.
		'''

		self.filename=filename
		self.readfile()

	def readfile(self):
		'''Force reread of the file and setting the data'''
		self.data=""
		try:
			with open (self.filename, "r") as myfile:
				self.data=myfile.read().replace('\n', '') #read the file and remove newline from the text
		except:
			pass # does nothing if error occurs


	def writefile(self, data, mode='w'):
		'''File may be updated by completely rewritting the file contents or appending the data to the end of the file.
		this is achieved by calling writefile(data, mode) method, where data is the string data to be written and
		mode can be 'a' for append and 'w' for writting the file anew.
		'''
		with open (self.filename, mode) as myfile:
			myfile.write(data)

	def getText(self):
		'''
		Method getText() or calling the object itself returns string of data
		'''
		return self.data

	def __str__(self):
		'''
		Method getText() or calling the object itself returns string of data
		'''
		return self.getText()

class Tape:
	'''
	Special class that manages configuration of trisurf (commonly named tape). It can read and parse configuration from disk or parse it from string.
	'''

	def __init__(self):
		'''The object is instatiated by calling Tape() constructor without parameters'''
		return

	def readTape(self, tape='tape'):
		'''
		Tape is read and parsed by calling the readTape() method with optional tape parameter, which is full path to filename where the configuration is stored.
		If the tape cannot be read, it prints error and exits.
		'''
		try:
			self.config=configobj.ConfigObj(tape)
		except:
			print("Error reading or parsing tape file!\n")
			exit(1)

	def setTape(self, string):
		'''Method setTape(string) parses the string in memory that hold the tape contents.'''
		self.config=configobj.ConfigObj(io.StringIO(string))
		return

	def getValue(self,key):
		'''Method getValue(key) returns value of a single parsed setting named "key".'''

		return self.config[key]

	def __str__(self):
		'''Calling the object itself, it recreates the tape contents from parsed values in form of key=value.'''
		retval=""
		for key,val in self.config.iteritems():
			retval=retval+str(key)+" = "+str(val)+"\n"
		return retval



class Directory:
	'''
	Class deals with the paths where the simulation is run and data is stored.
	'''
	def __init__(self, maindir=".", simdir=""):
		'''Initialization Directory() takes two optional parameters, namely maindir and simdir. Defaults to current directory. It sets local variables maindir and simdir accordingly.'''
		self.maindir=maindir
		self.simdir=simdir
		return

	def fullpath(self):
		'''
		Method returns string of path where the data is stored. It combines values of maindir and simdir as maindir/simdir on Unix.
		'''
		return os.path.join(self.maindir,self.simdir)

	def exists(self):
		''' Method checks whether the directory  specified by fullpath() exists. It return True/False on completion.'''
		path=self.fullpath()
		if(os.path.exists(path)):
			return True
		else:
			return False

	def make(self):
		''' Method make() creates directory. If it fails it exits the program with error message.'''
		try:
			os.makedirs(self.fullpath())
		except:
			print("Cannot make directory "+self.fullpath()+"\n")
			exit(1)
		return

	def makeifnotexist(self):
		'''Method makeifnotexist() creates directory if it does not exist.'''
		if(self.exists()==0):
			self.make()
		return

	def remove(self):
		'''Method remove() removes directory recursively. WARNING! No questions asked.'''
		if(self.exists()):
			try:
				os.rmdir(self.fullpath())
			except:
				print("Cannot remove directory "+self.fullpath()+ "\n")
				exit(1)
		return

	def goto(self):
		'''
		Method goto() moves current directory to the one specified by fullpath(). WARNING: when using the relative paths, do not call this function multiple times.
		'''
		try:
			os.chdir(self.fullpath())
		except:
			print("Cannot go to directory "+self.fullpath()+"\n")
		return


class Statistics:
	'''
	Class that deals with the statistics file from the simulations.
	File is generally large and not all data is needed, so it is dealt with in a specific way.
	'''

	def __init__(self,path,filename="statistics.csv"):
		'''
		At the initialization call it receives optional filename parameter specifying the path and filename of the statistics file.

		The local variables path, filename, fullname (joined path and filename) and private check if the file exists are stored.
		'''
		self.path=path
		self.filename=filename
		self.fullname=os.path.join(path,filename)
		self.fileOK=self.read()
		return

	def exists(self):
		'''Method check if the statistics file exists.'''
		if(os.path.isfile(self.fullname)):
			return True
		else:
			return False

	def mapcount(self):
		'''
		Internal method for determining the number of the lines in the most efficient way. Is it really the most efficient?
		'''
		f = open(self.fullname, "r+")
		buf = mmap.mmap(f.fileno(), 0)
		lines = 0
		readline = buf.readline
		while readline():
			lines += 1
		return lines

	def read(self):
		'''
		Method read() reads the statistics if it exists. It sets local variable dT storing the time differential between two intervals of simulation (outer loops). It also stores last simulation loop and the start of the run.
		'''
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
				#print("Cannot read statistics file in "+self.fullname+"\n")
				return(False)
		else:
			#print("File "+self.fullname+" does not exists.\n")
			return(False)

		self.dT=(int(epoch2)-int(epoch1))/(int(n2)-int(n1))
		self.last=n2
		self.startDate=epoch1
		return(True)

	def __str__(self):
		'''
		Prints the full path with filename of the statistics.csv file
		'''
		return(str(self.fullname))



class Runner:
	'''
	Class Runner consists of a single running or terminated instance of the trisurf. It manages starting, stopping, verifying the running process and printing the reports of the configured instances.
	'''
	def __init__(self, subdir='run0', tape='', snapshot=''):
		self.subdir=subdir
		self.fromSnapshot=False
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
		self.fromSnapshot=True
		self.snapshotFile=snapshotfile
		root = tree.getroot()
		tapetxt=root.find('tape')
		version=root.find('trisurfversion')
		self.tape=Tape()
		self.tape.setTape(tapetxt.text)

	def getPID(self):
		self.Dir=Directory(maindir=self.maindir,simdir=self.subdir)
		self.Dir.makeifnotexist()
		try:
			fp = open(os.path.join(self.Dir.fullpath(),'.lock'))
		except IOError as e:
			return 0 #file probably does not exist. e==2??
		pid=fp.readline()
		fp.close()
		return int(pid)

	def getStatus(self):
		pid=self.getPID()
		if(pid==0):
			return TS_NOLOCK
		if(psutil.pid_exists(int(pid))):
			proc= psutil.Process(int(pid))
			if proc.name=="trisurf":
				if proc.status=="stopped":
					return TS_STOPPED
				else:
					return TS_RUNNING
			else:
				return TS_NONEXISTANT
		else:
			return TS_NONEXISTANT

	def start(self):
		if(self.getStatus()==0):
			self.Dir=Directory(maindir=self.maindir,simdir=self.subdir)
			self.Dir.makeifnotexist()
			cwd=Directory(maindir=os.getcwd())
			self.Dir.goto()
			print("Starting trisurf-ng executable at "+self.Dir.fullpath()+"\n")
			if(self.fromSnapshot==True):
				params=["trisurf", "--restore-from-vtk",self.snapshotFile]
			else:
				params="trisurf"
			Popen (params, stdout=False)
			cwd.goto()
		else:
			print("Process already running. Not starting\n")
		return

	def stop(self):
		pass

	def setMaindir(self,prefix,variables):
		maindir=""
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
		self.Dir=Directory(maindir=self.maindir,simdir=self.subdir)
		self.statistics=Statistics(self.Dir.fullpath(), statfile)
		self.Comment=FileContent(os.path.join(self.Dir.fullpath(),".comment"))
		pid=self.getPID();
		status=self.getStatus()
		if(status==TS_NONEXISTANT or status==TS_NOLOCK):
			statustxt="Not running"
			pid=""
		elif status==TS_STOPPED:
			statustxt="Stopped"
		else:
			statustxt="Running"

		if(self.statistics.fileOK):
			report=[time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(self.statistics.startDate))),str(datetime.timedelta(microseconds=(int(self.tape.config['iterations'])-int(self.statistics.last))*self.statistics.dT)*1000), statustxt, pid, str(self.Dir.fullpath()), self.Comment.getText()]
		else:
			report=["N/A","N/A\t",statustxt, pid, str(self.Dir.fullpath()), self.Comment.getText()]
		return report

	def writeComment(self, data, mode='w'):
		self.Dir=Directory(maindir=self.maindir,simdir=self.subdir)
		self.Comment=FileContent(os.path.join(self.Dir.fullpath(),".comment"))
		self.Comment.writefile(data,mode=mode)

	def __str__(self):
		if(self.getStatus()==0):
			str=" not running."
		else:
			str=" running."
		return(self.Dir.fullpath()+str)


