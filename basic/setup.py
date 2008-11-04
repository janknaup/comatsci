#!/usr/bin/env python
from distutils.core import setup,Extension
import os

startpath=os.getcwd()

VERSIONTAG="1.0.0"
AUTHOR="Jan M. Knaup"
AU_EMAIL="Knaup@bccms.uni-bremen.de"
URL="http://www.bccms.uni-bremen.de/en/people/home/j_m_knaup/software/"

setup (	name="comatsci-base",
		version=VERSIONTAG,
		packages=['comatsci', 'comatsci.Geometry','geostatspack'],
		package_dir={'comatsci':'src/comatsci','geostatspack':'src/geostatspack'},
		ext_package='comatsci',
		ext_modules=[Extension('Geometry.geoext',['src/extensions/geoext.c']),
			Extension('splext',['src/extensions/splext.c'])],
		scripts=['src/scripts/geostats','src/scripts/geoconv',
                         'src/scripts/coordination_check','src/scripts/splresample',
                         'src/scripts/splderive',],
		description="Basic Computational Materials Science Toolkit",
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
