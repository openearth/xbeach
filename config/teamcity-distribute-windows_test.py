import os
	
def distribute(scptarget,revision,bindir,workdir):
	import glob
	import shutil
		
	releases = {}
	releases["Release"] = ["trunk version","trunk"]
	releases["netcdf_Release"] = ["with netcdf support","netcdf"]
	releases["MPI_Release"] = ["with mpi support","mpi"]
	releases["MPI_netcdf_Release"] = ["with mpi and netcdf support","mpi_netcdf"]

	pattern = '%s*/*/*.zip'%(bindir)
	files = glob.glob(pattern)

	ln = ''
	
	# file = files[0]
	#for file in files:
	#	nameparts = file.replace('\\','/').replace(bindir.replace('\\','/'),"").split('/')
	#	platform = nameparts[0]
	#	configuration = nameparts[1]
	#	filename = nameparts[2]
		
	#	# Copy file to bin dir (with correct name)
	#	newname = "xbeach_" + str(revision) + "_" + platform + "_" + releases[configuration][1] + ".zip"
	#	scpcopyfile(workdir, file, scptarget + "bin/" + newname)
		
	#	# fill template list string
	#	label = "XBeach rev. %d %s (%s)" % (revision,platform,releases[configuration][0])
	#	ln += '<li><a href="bin/%s">%s</a></li>%s' % (newname,label,os.linesep)
		
	# string replace in tmp file
	workhtml = workdir + "index_bin.html"
	#with open(workhtml, "wt") as fout:
	#	with open(workdir + "index_bin.html.tmpl", "rt") as fin:
	#		for line in fin:
	#		    fout.write(line.replace('${list}', ln))
	
	# Copy all zips and html to oss site
	scpcopyfile(workdir, workhtml, scptarget + "index_bin.html")

def scpcopyfile(workdir, source,destination):
	import subprocess
	import os
	
	FNULL = open(os.devnull,'w')
	pscp = workdir + "pscp.exe"
	
	args = pscp + " -l os.environ['XBeachOssUser'] -pw os.environ['XBeachOssPass'] -r -v " + source + " " + destination
	
	subprocess.call(args, stdout=FNULL, stderr=FNULL, shell=False)

revision = os.environ['RevisionNumber']
scptarget = "oss.deltares.nl:/drbd/www/content/xbeach/testbed/"
bindir = os.environ['XBeachDistributeDir']
workdir = os.environ['XBeachWorkingDir']

distribute(scptarget,revision,bindir,workdir)