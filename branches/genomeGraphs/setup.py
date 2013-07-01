from distutils.core import setup
import os

"""
_MEMO_ For a fresh re-install

cd /go//to/setup.py

python setup.py sdist
cd dist/
tar zxvf genomeGraphs-0.1.0.tar.gz
cd genomeGraphs-0.1.0
python setup.py install --user --install-scripts $HOME/bin/
cd ../../
"""

setup(
   name = 'genomeGraphs',
   version = '0.1.0',
   author = 'Dario Beraldi',
   author_email = 'dario.beraldi<at>gmail<dot>com',
   url= '',
   description = "Command-line oriented genome viewer",
   classifiers = [
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License (GPL)',
      'Operating System :: POSIX',
      'Programming Language :: Python'
   ],

   requires = [ 'pybedtools', 'python (>=2.6, <3.0)' ],
   
   py_modules = [
      'genomeGraphs.genomeGraphs',
      'genomeGraphs.pycoverage',
      'genomeGraphs.pympileup',
      'genomeGraphs.validate_args'
   ],

   scripts = [
      'scripts/genomeGraphs',
   ],
   packages= ['genomeGraphs'],
   package_data= {'genomeGraphs': ['R_template.R', 'mpileupToNucCounts.jar']}
#   data_files=[
#      ('genomeGraphs', ['genomeGraphs/R_template.R']),
#      ('genomeGraphs/java_code/mpileupToNucCounts',  ['genomeGraphs/java_code/mpileupToNucCounts/Pile.class',
#                                                      'genomeGraphs/java_code/mpileupToNucCounts/mpileupParser.class',
#                                                      'genomeGraphs/java_code/mpileupToNucCounts/mpileupParser.java'])
#   ]
)
