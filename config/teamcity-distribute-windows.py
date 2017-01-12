import os
import sys

def distribute(scptarget,revision,bindir,workdir,user,password):
	import glob
	import shutil
		
	releases = {}
	releases["Release"] = ["trunk version","trunk"]
	releases["netcdf_Release"] = ["with netcdf support","netcdf"]
	releases["MPI_Release"] = ["with mpi support","mpi"]
	releases["MPI_netcdf_Release"] = ["with mpi and netcdf support","mpi_netcdf"]

	pattern = '%s*/*/*.zip'%(bindir)
	files = glob.glob(pattern)
	
	print "##teamcity[message text='Copying binaries']"
	
	ln = ''
	for file in files:
		nameparts = file.replace('\\','/').replace(bindir.replace('\\','/'),"").split('/')
		platform = nameparts[0]
		configuration = nameparts[1]
		filename = nameparts[2]
		
		# Copy file to bin dir (with correct name)
		newname = "xbeach_" + str(revision) + "_" + platform + "_" + releases[configuration][1] + ".zip"
		
		scpcopyfile(workdir, file, scptarget + "bin/" + newname,user,password)
		
		# fill template list string
		label = "XBeach rev. %d %s (%s)" % (revision,platform,releases[configuration][0])
		ln += '<li><a href="bin/%s">%s</a></li>%s' % (newname,label,os.linesep)
		
	print "##teamcity[message text='Finished copying binaries']"
	
	print "##teamcity[message text='Creating index html of download page']"
	# string replace in tmp file
	workhtml = workdir + "\\index_bin.html"
	with open(workhtml, "wt") as fout:
		with open(workdir + "\\index_bin.html.tmpl", "rt") as fin:
			for line in fin:
			    fout.write(line.replace('${list}', ln))
	
	# Copy all zips and html to oss site
	scpcopyfile(workdir, workhtml, scptarget + "index_bin.html",user,password)

def scpcopyfile(workdir, source,destination,user,password):
	import subprocess
	import os
	
	FNULL = open(os.devnull,'w')
	pscp = workdir + "\\pscp.exe"
	
	args = "echo y | " + pscp + " -l " + user + " -pw " + password + " -r -v " + source + " " + destination
	
	print "##teamcity[message text='Copying %s']"%(source)
	subprocess.call(args, stdout=FNULL, stderr=FNULL, shell=True)

if not(len(sys.argv)==6):
	raise Exception('Wrong number of input arguments')
	
revision = int(sys.argv[1])
scptarget = "v-oss003.dlt.proteon.nl:/htdocs/xbeach/testbed/"
bindir = sys.argv[2]
workdir = sys.argv[3]
user = sys.argv[4]
password = sys.argv[5]

print " "
print "Publish XBeach executables to www.xbeach.org"
print "  "
print "Revision number   : %i"%(revision)
print "distributes dir   : %s"%(bindir)
print "Work directory    : %s"%(workdir)
print "  "

distribute(scptarget,revision,bindir,workdir,user,password)
