import sys, getopt
import tabulate

def start(Runs):
	argv=sys.argv[1:]
	processno=0
	try:
		opts, args = getopt.getopt(argv,"a:n:hrsc:")
	except getopt.GetoptError:
		print('tsmgr [-n process number] [-h] [-r] [-s] [-c comment text] [-a comment text]')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print ('tsmgr [-n process number] [-h] [-r] [-s] [-c comment text] [-a comment text]')
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
				#print(reportstr)
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
			print('tsmgr [-n process number] [-h] [-r] [-s] [-c comment text] [-a comment text]')
			sys.exit(2)






