import sys, getopt
import tabulate


def printHelp():
	print('Python module tsmgr accept following switches:\n')
	print('tsmgr [-n process number] [-R] [-h] [-r] [-s] [-c comment text] [-a comment text]\n')
	print('[-n process number]: number of process for which -s -r -c or -a switch apply. Should be placed before any other switch');
	print('[-R]               : raw output for -s switch');
	print('[-r]               : run process');
	print('[-s]               : process status');
	print('[-c comment text]  : write new comment for process');
	print('[-a comment text]  : append additional comment for process');
	print('[-h]               : print help');


def start(Runs):
	argv=sys.argv[1:]
	processno=0
	raw=False
	try:
		opts, args = getopt.getopt(argv,"Ra:n:hrsc:")
	except getopt.GetoptError:
		printHelp()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-R':
			raw=True
		elif opt == '-h':
			printHelp()
			sys.exit()
		elif opt == '-r':
			if processno:
				localRuns=[Runs[processno-1]]
			else:
				localRuns=Runs
			for run in localRuns:
				run.start()
		elif opt == '-s':
			report=[]
			i=1
			if processno:
				localRuns=[Runs[processno-1]]
			else:
				localRuns=Runs
			for run in localRuns:
				line=run.getStatistics()
				line.insert(0,i)
				report.append(line)
				i=i+1
			if(raw):
				print(report)
			else:
				print ("\n\nTrisurf running processes report\n")
				print (tabulate.tabulate(report,headers=["Run no.", "Run start time", "ETA", "Status", "PID", "Path", "Comment"], tablefmt='fancy_grid'))
		elif opt == '-n':
			processno=int(arg)
			if processno<1 or processno>len(Runs) :
				processno=0
		elif opt == '-c':
			comment = arg
			if processno:
				Runs[processno-1].writeComment(arg)
		elif opt == '-a':
			comment = arg
			if processno:
				Runs[processno-1].writeComment("\n"+arg, 'a')

			
		else:
			printHelp()
			sys.exit(2)






