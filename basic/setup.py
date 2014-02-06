#!/usr/bin/env python
from setuptools import setup,Extension
import os,sys

execfile('src/comatsci/_version.py')


VERSIONTAG=VERSION

startpath=os.getcwd()

AUTHOR="Jan M. Knaup"
AU_EMAIL="janknaup@gmail.com"
URL="http://www.bccms.uni-bremen.de/en/people/home/j_m_knaup/software/"

# define a list of script aliases, these will be installed either by copying
# the source script or by symlinking, if symlinks are available on the target
# platform, format is a dictionary with script names as keys containing lists of
# alias name strings
aliases={"comatsci.geoconvert":["geoconv","togen","toxyz","tofmg","tofdf","topdb","toxyzq","totm","toaims","tocdh","tovasp"]}

aliasentries=[]
for cmd in aliases.keys():
 		for alias in aliases[cmd]:
 			aliasentries.append("{0:s} = {1:s}:mainfunc".format(alias,cmd))

print aliases
print aliasentries

#sys.exit() 			

distrib=setup (	name="comatsci-base",
		version=VERSIONTAG,
		packages=['comatsci', 'comatsci.geometry',"comatsci.path","comatsci.optimizers",
			  "comatsci.schedulers","comatsci.calculators","comatsci.calculators.potentials"],
		package_dir={'comatsci':'src/comatsci',},
		ext_package='comatsci',
		ext_modules=[Extension('geometry.geoext',['src/extensions/geoext.c'],
							extra_compile_args=['-fopenmp','-funroll-loops','-O2'],extra_link_args=['-lgomp']),
			Extension('splext',['src/extensions/splext.c'])],
		scripts=['src/scripts/coordination_check','src/scripts/splresample',
                         'src/scripts/splderive','src/scripts/dumpbonds','src/scripts/fitrep',
			 "src/scripts/scale_linkdists","src/scripts/chargeanalys-2D",
			 "src/scripts/dosanalys-2D","src/scripts/dosanalys-3D",
			 'src/scripts/pastafarian','src/scripts/pathprepare',
			'src/scripts/pathprops','src/scripts/multiaverage',
			'src/scripts/dosplot','src/scripts/pathconvert'],
		data_files=[("share/doc/comatsci",["doc/comatsci.pdf",
						   "doc/run.sh.example"])],
		description="Computational Materials Science Toolkit",
		author=AUTHOR,
		author_email=AU_EMAIL, 
		url=URL, 
		install_requires=["h5py >= 2.0","numpy >= 1.6",],
		classifiers=[
			'Development Status :: 5 - Production/Stable',
			'Environment :: Console',
			'Intended Audience :: Education',
			'Intended Audience :: Developers',
			'Intended Audience :: Science/Research',
			'License :: OSI Approved :: Non-Profit Open Software License v3',
			'Operating System :: POSIX',
			'Programming Language :: Python',
			'Programming Language :: Python :: 2.7',
			'Topic :: Scientific/Enginieering',
			'Topic :: Scientific/Enginieering :: Physics',
			'Topic :: Scientific/Enginieering :: Chemistry',
		],
		entry_points = {
        'console_scripts': []+aliasentries,
        'gui_scripts': []
    	}
		)

# #********************************************************************************************
# # post-install operations
# #********************************************************************************************
# 
# # scripts postprocessing - check if scripts were installed
# 
# if distrib.have_run.get("install_scripts",0)==1:
# 	# script aliases processing
# 	# decide wheter to symlink or copy aliases
# 	if "symlink" in dir(os):
# 		linktype="sym"
# 	else:
# 		linktype=None
# 	# store install target and copy command locally to save typing
# 	installdir=os.path.abspath(distrib.command_obj["install_scripts"].install_dir)
# 	# process aliases
# 	os.chdir(installdir)
# 	for cmd in aliases.keys():
# 		for alias in aliases[cmd]:
# 			# remove existing alias if present
# 			if os.path.exists(alias) and not distrib.dry_run:
# 				os.unlink(alias)
# 			# create new alias
# 			copycmd=distrib.command_obj["install_scripts"].copy_file(cmd,alias,link=linktype)
# 	os.chdir(startpath)
# 	# finished with script aliases
# # finished with scripts postprocessing
# 
# #print distrib.command_obj
