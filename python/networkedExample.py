from trisurf import trisurf
from trisurf import tsmgr




#Ok... Configure your keys:
#ssh-keygen
#and copy them to all the remote hosts
#ssh-copy-id -i ./ssh/id_rsa.pub username@remotehost

run2=trisurf.Runner(tape='tape')
run2.setMaindir(("N","k","V","Np","Nm"),("nshell","xk0","constvolswitch","npoly","nmono"))
run2.setSubdir("run1")

run3=trisurf.Runner(tape='tape')
run3.setMaindir(("N","k","V","Np","Nm"),("nshell","xk0","constvolswitch","npoly","nmono"))
run3.setSubdir("run2")



Runs=[run2, run3]

hosts=({'name':'natalie','address':'kabinet.penic.eu', 'runs':Runs, 'username':'samo'},
	{'name':'altea','address':'127.0.0.1', 'runs':Runs, 'username':'samo'})

tsmgr.start(hosts)
