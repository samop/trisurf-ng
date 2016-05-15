#!/usr/bin/python3
from trisurf import tsmgr
from trisurf import trisurf



print("Running trisurf version "+ tsmgr.getTrisurfVersion())

#Simple example how to start simulation from a previos snapshot
run1=trisurf.Runner(snapshot='snapshot.vtu')
run1.setMaindir(("N","k","V","Np","Nm"),("nshell","xk0","constvolswitch","npoly","nmono"))
run1.setSubdir("run0")

#Example how to start simulation from tape. Extra argument in runArgs will be passed to trisurf executable (meaning that simulation will always start from the beginning (bipyramid) ignoring the fact that some states may have been calculated already)
run2=trisurf.Runner(tape='tape', runArgs=['--force-from-tape'])
run2.setMaindir(("N","k","V","Np","Nm"),("nshell","xk0","constvolswitch","npoly","nmono"))
run2.setSubdir("run1")

#Example of programatical setup of 4 runs
pRun=[]
for i in range(0,4): #0,1,2,3
	tpRun=trisurf.Runner(tape='tape')
	tpRun.setMaindir(("N","k","V","Np","Nm"),("nshell","xk0","constvolswitch","npoly","nmono"))
	tpRun.setSubdir("programatical_"+str(i))
	pRun.append(tpRun)


#obligatory final configuration step: combine all runs
Runs=[run1,run2]+pRun
#start manager with configured runs
tsmgr.start(Runs)

