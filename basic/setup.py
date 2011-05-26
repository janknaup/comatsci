#!/usr/bin/env python
from distutils.core import setup,Extension
from distutils.dir_util import remove_tree
import os

startpath=os.getcwd()

VERSIONTAG="1.1.0-rc2"
AUTHOR="Jan M. Knaup"
AU_EMAIL="Knaup@bccms.uni-bremen.de"
URL="http://www.bccms.uni-bremen.de/en/people/home/j_m_knaup/software/"

# define a list of script aliases, these will be installed either by copying
# the source script or by symlinking, if symlinks are available on the target
# platform, format is a dictionary with script names as keys containing lists of
# alias name strings
aliases={"geoconv":["togen","toxyz","tofmg","tofdf","topdb","toxyzq","totm"]}

distrib=setup (	name="comatsci-base",
		version=VERSIONTAG,
		packages=['comatsci', 'comatsci.Geometry'],
		package_dir={'comatsci':'src/comatsci',},
		ext_package='comatsci',
		ext_modules=[Extension('Geometry.geoext',['src/extensions/geoext.c']),
			Extension('splext',['src/extensions/splext.c'])],
		scripts=['src/scripts/geoconv',
                         'src/scripts/coordination_check','src/scripts/splresample',
                         'src/scripts/splderive','src/scripts/dumpbonds'],
		data_files=[("share/doc/comatsci",["doc/comatsci-basic.pdf"])],
		description="Basic Computational Materials Science Toolkit",
		author=AUTHOR,
		author_email=AU_EMAIL, 
		url=URL, 
		classifiers=[
			'Development Status :: 5 - Production/Stable',
			'Environment :: Console',
			'Intended Audience :: Education',
			'Intended Audience :: Developers',
			'Intended Audience :: Science/Research',
			'License :: OSI Approved :: Non-Profit Open Software License v3',
			'Operating System :: POSIX',
			'Programming Language :: Python',
			'Programming Language :: Python :: 2.4',
			'Programming Language :: Python :: 2.5',
			'Topic :: Scientific/Enginieering',
			'Topic :: Scientific/Enginieering :: Physics',
			'Topic :: Scientific/Enginieering :: Chemistry',
		]
		)

#********************************************************************************************
# post-install operations
#********************************************************************************************

# scripts postprocessing - check if scripts were installed

if distrib.have_run.get("install_scripts",0)==1:
	# script aliases processing
	# decide wheter to symlink or copy aliases
	if "symlink" in dir(os):
		linktype="sym"
	else:
		linktype=None
	# store install target and copy command locally to save typing
	installdir=os.path.abspath(distrib.command_obj["install_scripts"].install_dir)
	# process aliases
	os.chdir(installdir)
	for cmd in aliases.keys():
		for alias in aliases[cmd]:
			# remove existing alias if present
			if os.path.exists(alias) and not distrib.dry_run:
				os.unlink(alias)
			# create new alias
			copycmd=distrib.command_obj["install_scripts"].copy_file(cmd,alias,link=linktype)
	os.chdir(startpath)
	# finished with script aliases
# finished with scripts postprocessing

#print distrib.command_obj
