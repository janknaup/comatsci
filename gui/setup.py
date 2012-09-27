#!/usr/bin/env python
from distutils.core import setup,Extension
from distutils.dir_util import remove_tree
import os,sys

sys.path.append("../basic/src/comatsci")

from constants import VERSIONPREFIX as COMATSCI_VERSIONPREFIX
from constants import VERSION as COMATSCI_VERSION

try:
	from bzrlib.branch import BzrBranch
	branch = BzrBranch.open_containing('.')[0]
	if branch.nick!="trunk":
		VERSIONTAG=COMATSCI_VERSIONPREFIX+"-{1:s}-{0:d}".format(branch.last_revision_info()[0],branch.nick)
	else:
		VERSIONTAG=COMATSCI_VERSIONPREFIX+"-{0:d}".format(branch.last_revision_info()[0])
except:
		VERSIONTAG=COMATSCI_VERSIONPREFIX

startpath=os.getcwd()


AUTHOR="Jan M. Knaup"
AU_EMAIL="janknaup@gmail.com"
URL="http://www.bccms.uni-bremen.de/en/people/home/j_m_knaup/software/"

# define a list of script aliases, these will be installed either by copying
# the source script or by symlinking, if symlinks are available on the target
# platform, format is a dictionary with script names as keys containing lists of
# alias name strings

distrib=setup (	name="comatsci-gui",
		version=VERSIONTAG,
		packages=['geostatspack4'],
		package_dir={'geostatspack4':'src/geostatspack4'},
		scripts=['src/scripts/geostats',],
		data_files=[("share/doc/comatsci",["doc/comatsci-gui.pdf"])],
		description="Computational Materials Science Toolkit GUI utilities",
		author=AUTHOR,
		author_email=AU_EMAIL, 
		url=URL, 
		classifiers=[
			'Development Status :: 5 - Production/Stable',
			'Environment :: X11 Applications :: Qt',
			'Intended Audience :: Education',
			'Intended Audience :: Developers',
			'Intended Audience :: Science/Research',
			'License :: OSI Approved :: Non-Profit Open Software License v3',
			'Operating System :: POSIX',
			'Programming Language :: Python',
			'Programming Language :: Python :: 2.6',
			'Programming Language :: Python :: 2.7',
			'Topic :: Scientific/Enginieering',
			'Topic :: Scientific/Enginieering :: Physics',
			'Topic :: Scientific/Enginieering :: Chemistry',
		]
		)

#********************************************************************************************
# post-install operations
#********************************************************************************************

#print distrib.command_obj
