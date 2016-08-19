import sys, platform
import os, shutil, hashlib, re, fnmatch

"""Description

Setup script for CPAT  -- Coding Potential Assessment Tool

Copyright (c) 2012 Liguo Wang <wangliguo78@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).
"""

if sys.version_info[0] != 2 or sys.version_info[1] < 7:
	print >> sys.stderr, "ERROR: CPAT only works on Python2.x and requires Python 2.7 or greater"
	sys.exit()

# Automatically download setuptools if not available
from distribute_setup import use_setuptools
use_setuptools()

from setuptools import *
from glob import glob

try:
	import numpy
	have_numpy = True
except:
	have_numpy = False
       
def main():
	setup(  name = "CPAT",
            version = "1.2.2",
            py_modules = [ 'psyco_full' ],
            packages = find_packages( 'lib' ),
            package_dir = { '': 'lib' },
            package_data = { '': ['*.ps'] },
            scripts = glob( "bin/*.py"),
            ext_modules = get_extension_modules(),
            test_suite = 'nose.collector',
            setup_requires = ['nose>=0.10.4','cython>=0.12'],
            author = "Liguo Wang",
            author_email ="wangliguo78@gmail.com",
	    platforms = ['Linux','MacOS'],
            description = "Coding Potential Assessment Tool",
            url = "http://rna-cpat.sourceforge.net/",
            zip_safe = False,
            dependency_links = ['http://bitbucket.org/james_taylor/bx-python/wiki/Home','http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/contents.html'],
            cmdclass=command_classes )

# ---- Commands -------------------------------------------------------------

from distutils.core import Command

# Use build_ext from Cython
command_classes = {}

# Use build_ext from Cython if found
try:
	import Cython.Distutils
	command_classes['build_ext'] = Cython.Distutils.build_ext
except:
	pass
    
# Use epydoc if found
try:
	import pkg_resources
	pkg_resources.require( "epydoc" )
	import epydoc.cli, sys, os, os.path
	# Create command class to build API documentation
	class BuildAPIDocs( Command ):
		user_options = []
		def initialize_options( self ):
			pass
		def finalize_options( self ):
			pass
		def run( self ):
			# Save working directory and args
			old_argv = sys.argv
			old_cwd = os.getcwd()
			# Build command line for Epydoc
			sys.argv = """epydoc.py bx --verbose --html --simple-term
                                       --exclude=._
                                       --exclude=_tests
                                       --docformat=reStructuredText
                                       --output=../doc/docbuild/html/apidoc""".split()
			# Make output directory
			if not os.path.exists( "./doc/docbuild/html/apidoc" ):
				os.mkdir( "./doc/docbuild/html/apidoc" )
			# Move to lib directory (so bx package is in current directory)
			os.chdir( "./lib" )
			# Invoke epydoc
			epydoc.cli.cli()
			# Restore args and working directory
			sys.argv = old_argv
			os.chdir( old_cwd )
	# Add to extra_commands    
	command_classes['build_apidocs'] = BuildAPIDocs
except:
	pass

# ---- Extension Modules ----------------------------------------------------

def get_extension_modules():
	extensions = []
	# Bitsets
	extensions.append( Extension( "bx.bitset",
                                  [ "lib/bx/bitset.pyx", 
                                    "src/binBits.c",
                                    "src/kent/bits.c",
                                    "src/kent/common.c" ],
                                  include_dirs=[ "src/kent", "src"] ) )
	# Interval intersection
	extensions.append( Extension( "bx.intervals.intersection", [ "lib/bx/intervals/intersection.pyx" ] ) )
	# Alignment object speedups
	extensions.append( Extension( "bx.align._core", [ "lib/bx/align/_core.pyx" ] ) )
	# NIB reading speedups
	extensions.append( Extension( "bx.seq._nib", [ "lib/bx/seq/_nib.pyx" ] ) )
	# 2bit reading speedups
	extensions.append( Extension( "bx.seq._twobit", [ "lib/bx/seq/_twobit.pyx" ] ) )
	# Translation if character / integer strings 
	extensions.append( Extension( "bx._seqmapping", [ "lib/bx/_seqmapping.pyx" ] ) )
	# BGZF
	extensions.append( Extension( "bx.misc.bgzf",
                                  [ "lib/bx/misc/bgzf.pyx", "src/samtools/bgzf.c" ],
                                  include_dirs=[ "src/samtools"],
                                  libraries=['z'] ) )

    
	# The following extensions won't (currently) compile on windows
	if platform.system() not in ( 'Microsoft', 'Windows' ):
        
		# Interval clustering                
		extensions.append( Extension( "bx.intervals.cluster",
                                  [ "lib/bx/intervals/cluster.pyx", 
                                    "src/cluster.c"],
                                  include_dirs=["src"] ) )
		# Position weight matrices
		extensions.append( Extension( "bx.pwm._position_weight_matrix",
                                  [ "lib/bx/pwm/_position_weight_matrix.pyx", "src/pwm_utils.c" ],
                                  include_dirs=["src"]  ) )
 	
		if have_numpy:
			extensions.append( Extension( "bx.motif._pwm", [ "lib/bx/motif/_pwm.pyx" ], 
                                          include_dirs=[numpy.get_include()] ) )
            
			# Sparse arrays with summaries organized as trees on disk
			extensions.append( Extension( "bx.arrays.array_tree", [ "lib/bx/arrays/array_tree.pyx" ], include_dirs=[numpy.get_include()] ) )  
        
			# Reading UCSC "big binary index" files
			extensions.append( Extension( "bx.bbi.bpt_file", [ "lib/bx/bbi/bpt_file.pyx" ] ) )
			extensions.append( Extension( "bx.bbi.cirtree_file", [ "lib/bx/bbi/cirtree_file.pyx" ] ) )
			extensions.append( Extension( "bx.bbi.bbi_file", [ "lib/bx/bbi/bbi_file.pyx" ], include_dirs=[numpy.get_include()] ) )
			extensions.append( Extension( "bx.bbi.bigwig_file", [ "lib/bx/bbi/bigwig_file.pyx" ], include_dirs=[numpy.get_include()] ) )
			extensions.append( Extension( "bx.bbi.bigbed_file", [ "lib/bx/bbi/bigbed_file.pyx" ], include_dirs=[numpy.get_include()] ) )
	
		# Reading UCSC bed and wiggle formats
		extensions.append( Extension( "bx.arrays.bed", [ "lib/bx/arrays/bed.pyx" ] ) )
		extensions.append( Extension( "bx.arrays.wiggle", [ "lib/bx/arrays/wiggle.pyx" ] ) )
	
	
		# CpG masking
		extensions.append( Extension( "bx.align.sitemask._cpg", \
                                      [ "lib/bx/align/sitemask/_cpg.pyx", 
                                        "lib/bx/align/sitemask/find_cpg.c" ] ) )
        
		# Counting n-grams in integer strings
		extensions.append( Extension( "bx.intseq.ngramcount", [ "lib/bx/intseq/ngramcount.pyx" ] ) )
	
		# Seekable access to bzip2 files
		extensions.append( Extension( "bx.misc._seekbzip2", 
                                      [ "lib/bx/misc/_seekbzip2.pyx",
                                        "src/bunzip/micro-bunzip.c" ],
                                      include_dirs=[ "src/bunzip" ] ) )
	#samtools
	if platform.system()=='Windows':
		include_os = ['win32']
		os_c_files = ['win32/getopt.c']
	else:
		include_os = []
		os_c_files = []
	extensions.append(Extension("csamtools",                   # name of extension
                                      [ "lib/pysam/csamtools.pyx" ]  +\
                                      [ "lib/pysam/%s" % x for x in (
                                      "pysam_util.c", )] +\
                                      glob( os.path.join( "lib/samtools", "*.pysam.c" )) +\
                                      os_c_files + \
                                      glob( os.path.join( "lib/samtools", "*", "*.pysam.c" ) ),
                                      library_dirs=[],
                                      include_dirs=[ "lib/samtools", "lib/pysam" ] + include_os,
                                      libraries=[ "z", ],
                                      language="c",
                                      define_macros = [('_FILE_OFFSET_BITS','64'),
                                      ('_USE_KNETFILE','')], 
                                      ))

	extensions.append(Extension("ctabix",                   # name of extension
                                      [ "lib/pysam/ctabix.pyx", ]  +\
                                      [ "lib/pysam/%s" % x for x in ( "tabix_util.c", )] +\
                                      os_c_files + \
                                      glob( os.path.join( "lib/tabix", "*.pysam.c" ) ),
                                      library_dirs=[],
                                      include_dirs=[ "lib/tabix", "lib/pysam" ] + include_os,
                                      libraries=[ "z", ],
                                      language="c",
                                      ))
	extensions.append(Extension("TabProxies",                   # name of extension
                                      [ "lib/pysam/TabProxies.pyx", ] + os_c_files,
                                      library_dirs=[],
                                      include_dirs= include_os,
                                      libraries=[ "z", ],
                                      language="c",
                                      ))

	extensions.append(Extension("cvcf",                   # name of extension
                                      [ "lib/pysam/cvcf.pyx", ] + os_c_files,
                                      library_dirs=[],
                                      include_dirs= ["lib/tabix",] + include_os,
                                      libraries=[ "z", ],
                                      language="c",
                                      ))


	return extensions     
 
# ---- Monkey patches -------------------------------------------------------

def monkey_patch_doctest():
    #
    # Doctest and coverage don't get along, so we need to create
    # a monkeypatch that will replace the part of doctest that
    # interferes with coverage reports.
    #
    # The monkeypatch is based on this zope patch:
    # http://svn.zope.org/Zope3/trunk/src/zope/testing/doctest.py?rev=28679&r1=28703&r2=28705
    #
	try:
		import doctest
		_orp = doctest._OutputRedirectingPdb
		class NoseOutputRedirectingPdb(_orp):
			def __init__(self, out):
				self.__debugger_used = False
				_orp.__init__(self, out)

			def set_trace(self):
				self.__debugger_used = True
				_orp.set_trace(self)

			def set_continue(self):
				# Calling set_continue unconditionally would break unit test coverage
				# reporting, as Bdb.set_continue calls sys.settrace(None).
				if self.__debugger_used:
					_orp.set_continue(self)
		doctest._OutputRedirectingPdb = NoseOutputRedirectingPdb
	except:
		pass

def monkey_patch_numpy():
	# Numpy pushes its tests into every importers namespace, yeccch.
	try:
		import numpy
		numpy.test = None
	except:
		pass


if __name__ == "__main__":
	monkey_patch_doctest()
	if have_numpy:
		monkey_patch_numpy()
	main()
