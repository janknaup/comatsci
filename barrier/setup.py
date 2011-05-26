#!/usr/bin/env python
from distutils.core import setup,Extension
import os,sys

startpath=os.getcwd()

VERSIONTAG="1.1.0-rc1"
AUTHOR="Jan M. Knaup"
AU_EMAIL="Knaup@bccms.uni-bremen.de"
URL="http://www.bccms.uni-bremen.de/en/people/home/j_m_knaup/software/"

distrib=setup (	name="comatsci-barriers",
		version=VERSIONTAG,
		packages=['comatsci.Calculators','comatsci.Calculators.Potentials','comatsci.Path','comatsci.Schedulers'],
		package_dir={'comatsci':'src/comatsci',},
		scripts=['src/scripts/pastafarian','src/scripts/pathprepare',
			'src/scripts/pathprops','src/scripts/multiaverage'],
		data_files=[("share/doc/comatsci",["doc/comatsci-barrier.pdf"])],
		description="Computational Materials Science Toolkit - Reaction Barrier Search Module",
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

