import paramiko



class Connection:
	def __init__(self, hostname, port=22, username=None, password=None):
		self.hostname=hostname
		self.port=port
		if(username!=None):
			self.username=username
		else:
			self.username=''
		if(password!=None):
			self.password=password
		else:
			self.password=''
		self.ssh=paramiko.SSHClient()
		self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		self.connected=False
		return
	
	def connect(self, Timeout=5):
		if(not self.connected):
			try:
				print("Trying to connect to: "+self.username+"@"+self.hostname+":"+str(self.port)+".")
				self.ssh.connect(self.hostname, username=self.username, password=self.password, port=self.port, timeout=Timeout)
				self.connected=True
			except:
				print("Error establishing connection with "+self.username+"@"+self.hostname+":"+str(self.port)+".")
				exit(1)
		else:
			print("Already connected!")
		return

	def disconnect(self):
		if(self.connected):
			try:
				self.ssh.close()
			except:
				print("Cannot disconect. Unknown error.")
		else:
			print("Cannot disconect. Already disconnected.")
		self.connected=False

	def execute(self,command):
		if(self.connected):
			try:
				stdin,stdout,stderr=self.ssh.exec_command(command)
				output=stdout.readlines()
				errors=stderr.readlines()
			#	print(errors)
				return(output)
			except:
				print("Cannot execute remote commands")
		else:
			print("Cannot execute remote commands. Connect first.")

	def send_file(self, local, remote):
		sftp=self.ssh.open_sftp()
		sftp.put(local,remote)
		sftp.close()

	def receive_file(self,remote,local):
		sftp=self.ssh.open_sftp()
		sftp.get(remote,local)
		sftp.close()

	def mkdir_remote(self,directory):
		sftp=self.ssh.open_sftp()
		sftp.mkdir(directory)
		sftp.close()

	
