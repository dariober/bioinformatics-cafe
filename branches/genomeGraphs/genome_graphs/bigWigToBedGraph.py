#!/usr/bin/env python
import sys
import pybedtools
import os
# import bigWigToBedGraph
import platform
import tempfile
import subprocess
from distutils import spawn
#"""
#for exe in linux.x86_64 macOSX.i386 macOSX.ppc macOSX.x86_64
#do
#
#    wget http://hgdownload.cse.ucsc.edu/admin/exe/${exe}/bigWigToBedGraph
#    mv bigWigToBedGraph ${exe}.bigWigToBedGraph
#
#done
#"""

class BigWigToBedGraphException(Exception):
    pass

def getVersionForBw2bdg():

    ## First try to see if bigWigToBedGrpah is on PATH
    bwToBd= spawn.find_executable('bigWigToBedGraph')
    if bwToBd is not None:
        return(bwToBd)

    ## Try to use one of the exec:
    platf= platform.platform()
    proc= platform.processor()
    if 'linux' in platf.lower():
        myplat= 'linux'
    elif 'darwin' in platf.lower() or 'macos' in platf.lower():
        myplat= 'macOSX'
    else:
        raise BigWigToBedGraphException("Can't get platform")
    if 'i386' in proc:
        myproc= 'i386'
    elif 'x86_64' in proc:
        myproc= 'x86_64'
    elif 'ppc' in proc:
        myproc= 'ppc'
    else:
        raise BigWigToBedGraphException("Can't get processor")
    exe= myplat + '.' + myproc + '.' + 'bigWigToBedGraph'
    exe= os.path.abspath(os.path.join(os.path.split(bigWigToBedGraph.__file__)[0], exe))
    if not os.path.isfile(exe):
        raise BigWigToBedGraphException("Can't find bigWigToBedGraph")
    return(exe)

def bigWigToBedGraphExe(exe, inBigWig, inbed, tmpdir):
    """Convert bigWig to bedGraph at specified regions
    exe:
        Name and path to executable file bigWigToBedGraph.
    inBigWig:
        BW file to convert
    inbed:
        pybedTool obejct with regions to extract from bigWig
    tmpdir:
        Working dir. Inside this dir a tempdir will be created with inside the
        converted file.
    Return:
        Name and path to converted bedGraph file
    """
    outdir= tempfile.mkdtemp(prefix= 'bigWigToBedGraph_', dir= tmpdir)
    outBdg= os.path.join(outdir, os.path.split(inBigWig)[1])
    tmpbdg= os.path.join(tmpdir, 'tmp_bdg')
    minbed= inbed.sort().merge()
    for line in minbed:
        cmd= ''.join([exe, ' -chrom=', line.chrom, ' -start=', str(line.start), ' -end=', str(line.end), ' ', inBigWig, ' ', tmpbdg, '; cat ', tmpbdg, ' >> ', outBdg, '; rm ', tmpbdg])
        p= subprocess.Popen(cmd, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)
        stdout, stderr= p.communicate()
        if not p.returncode == 0:
            print(cmd)
            raise BigWigToBedGraphException(p.stderr)
    return(outBdg)

#inbed= pybedtools.BedTool('../actb.bed')
#outBdg= bigWigToBedGraphExe(getVersionForBw2bdg(), 'profile.odd.bw', inbed, tmpdir= '.')
#print(outBdg)
sys.exit()