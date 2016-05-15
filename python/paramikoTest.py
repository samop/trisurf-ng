#!/usr/bin/python3
import paramiko

ssh = paramiko.SSHClient()
#ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect('hostname', username='', password='',port=22,timeout=5)
stdin, stdout, stderr = ssh.exec_command('df -h')
x=stdout.readlines()
sftp = ssh.open_sftp()
#sftp.put(localpath, remotepath)
sftp.close()
ssh.close()
print (x)
