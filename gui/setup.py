#!/usr/bin/env python
from setuptools import setup,Extension
import os,sys


VERSIONTAG="1.6.0"

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
		data_files=[("share/doc/comatsci",["doc/comatsci-gui.pdf"])],
		description="Computational Materials Science Toolkit GUI utilities",
		author=AUTHOR,
		author_email=AU_EMAIL, 
		url=URL, 
		install_requires=["comatsci-base >= 1.4","matplotlib"],
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
		],
        entry_points = {
            'console_scripts': [],
            'gui_scripts': ["geostats = geostatspack4.geostats:mainfunc",]
    	}
)

