import os
import re
import sys
import subprocess
import paramiko
import zipfile


PLATFORMS = ['win32', 'x64']
CONFIGURATIONS = ['Release', 'netcdf_Release', 'MPI_Release', 'MPI_netcdf_Release']


class SSHClient:
    '''Dummy object to bypass paramiko'''

    def __init__(self, workdir):
        self.workdir = workdir
        
    
    def set_missing_host_key_policy(self, policy):
        self.missing_host_key_policy = policy
        
        
    def connect(self, host, port=22, username=None, password=None):
        self.host = host
        self.port = port
        self.username = username
        self.password = password
        
        
    def open_sftp(self):
        return self
        
        
    def put(self, src, dst):

        FNULL = open(os.devnull, 'w')
        pscp = self.workdir + '\\pscp.exe'
    
        args = 'echo y | %s -l %s -pw %s -r -v %s %s:%s' % (pscp, self.username, self.password, src, self.host, dst)
        subprocess.call(args, stdout=FNULL, stderr=FNULL, shell=True)


    def close(self):
        pass
    
    
def distribute(targetdir, revision, bindir, workdir, username, password, host='v-oss003.dlt.proteon.nl', port=22):
    '''Zip releases and upload to OSS website'''

    # open ssh socket
    ssh = SSHClient(workdir) #ssh = paramiko.SSHClient() # doesn't work because of buggy instalation on build agents
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host, port=port, username=username, password=password)
    sftp = ssh.open_sftp()

    print '##teamcity[message text=\'Uploading binaries\']'
    
    # loop over configurations
    ln = ''
    for platform in PLATFORMS:
        for configuration in CONFIGURATIONS:

            # create identifiers
            cid, clabel = get_configuration_identifiers(configuration)

            cpath = os.path.join(bindir, platform, configuration)
            if os.path.exists(cpath):
                
                # zip configuration files
                zname = 'xbeach_%d_%s_%s.zip' % (revision, platform, cid)
                print '##teamcity[message text=\'Creating archive %s\']' % zname
                zpath = zip_directory_contents(cpath, zname)

                # upload zip file
                print '##teamcity[message text=\'Uploading %s\']' % zname
                sftp.put(zpath, targetdir + 'bin/')
                
                # fill template list string
                label = 'XBeach rev. %d %s (%s)' % (revision, platform, clabel)
                ln += '<li><a href="bin/%s">%s</a></li>%s' % (zname, label, os.linesep)
        
    print '##teamcity[message text=\'Finished uploading binaries\']'
    
    print '##teamcity[message text=\'Creating index html of download page\']'

    # string replace in tmpl file (should be mako, but not available on build agents)
    htmlfile = os.path.join(workdir, 'index_bin.html')
    tmplfile = os.path.join(workdir, 'index_bin.html.tmpl')
    with open(htmlfile, 'wt') as fout:
        with open(tmplfile, 'rt') as fin:
            for line in fin:
                fout.write(line.replace('${list}', ln))
            
    # copy html to oss site
    sftp.put(htmlfile, targetdir)

    sftp.close()
    ssh.close()

    
def zip_directory_contents(cpath, zname):
    '''Zip entire directory and return path to resulting zip file

    Parameters
    ----------
    cpath : str
        Path to be zipped
    zname : str
        Name of zip file

    Returns
    -------
    zpath : str
        Path to zip file

    '''
    
    zpath = os.path.join(cpath, zname)
    with zipfile.ZipFile(zpath, 'w') as fp:
        for root, dirs, files in os.walk(cpath):
            for fname in files:
                if fname.endswith('.zip'):
                    continue
                fp.write(os.path.join(root, fname), fname)

    return zpath
                                                

def get_configuration_identifiers(configuration):
    '''Returns ID and label identifying a configuration'''

    cid = re.sub('_?Release$', '', configuration).lower()

    if not cid:
        cid = 'trunk'
        clabel = 'trunk version'
    else:
        clabel = 'with %s support' % ' and '.join(cid.split('_'))

    return cid, clabel


if __name__ == '__main__':

    if not(len(sys.argv)==6):
        raise Exception('Wrong number of input arguments')
    
    revision = int(sys.argv[1])
    targetdir = '/htdocs/xbeach/testbed/'
    bindir = sys.argv[2]
    workdir = sys.argv[3]
    username = sys.argv[4]
    password = sys.argv[5]

    print ' '
    print 'Publish XBeach executables to www.xbeach.org'
    print ' '
    print 'Revision number   : %i' % (revision)
    print 'Distribute dir    : %s' % (bindir)
    print 'Work directory    : %s' % (workdir)
    print ' '

    distribute(targetdir,
               revision,
               bindir,
               workdir,
               username,
               password)
