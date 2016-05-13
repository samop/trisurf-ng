from trisurf import tsmgr
from trisurf import trisurf

run1=trisurf.Runner(snapshot='snapshot.vtu')
run1.setMaindir(("N","k","V","Np","Nm"),("nshell","xk0","constvolswitch","npoly","nmono"))
run1.setSubdir("run0")

run2=trisurf.Runner(tape='tape', runArgs=['--force-from-tape'])
run2.setMaindir(("N","k","V","Np","Nm"),("nshell","xk0","constvolswitch","npoly","nmono"))
run2.setSubdir("run1")


#obligatory: combine all runs
Runs=[run1,run2]

tsmgr.start(Runs)

