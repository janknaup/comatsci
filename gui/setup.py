#!/usr/bin/env python
from distutils.core import setup,Extension
from distutils.dir_util import remove_tree
import os

startpath=os.getcwd()

VERSIONTAG="1.1.0-rc1"
AUTHOR="Jan M. Knaup"
AU_EMAIL="Knaup@bccms.uni-bremen.de"
URL="http://www.bccms.uni-bremen.de/en/people/home/j_m_knaup/software/"

# define a list of script aliases, these will be installed either by copying
# the source script or by symlinking, if symlinks are available on the target
# platform, format is a dictionary with script names as keys containing lists of
# alias name strings

distrib=setup (	name="comatsci-gui",
		version=VERSIONTAG,
		packages=['geostatspack'],
		package_dir={'geostatspack':'src/geostatspack'},
		scripts=['src/scripts/geostats',],
		data_files=[("share/doc/comatsci",["doc/comatsci-gui.pdf"])],
		description="Computational Materials Science Toolkit GUI utilities",
		author=AUTHOR,
		author_email=AU_EMAIL, 
		url=URL, 
		classifiers=[
			'Development Status :: 5 - Production/Stable',
			'Environment :: Console',
			'Environment :: X11 Applications :: Qt',
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

#print distrib.command_obj
