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
   version = '0.2.1a', ## MAKE IT MATCH genomeGraphs.py
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
      'genome_graphs.genomeGraphs',
      'genome_graphs.pycoverage',
      'genome_graphs.pympileup',
      'genome_graphs.validate_args'
   ],

   scripts = [
      'scripts/genomeGraphs',
      'scripts/genomeGraphsDemo.py',
   ],
    packages=['genome_graphs'],
    package_data = {
        'genome_graphs': ['R_template.R', 'mpileupToNucCounts.jar', 'demo/*'],
    },
)
