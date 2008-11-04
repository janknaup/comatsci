##############################################################################
# comatsci/__init__.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

#import constants definitions first
from constants import *

#import helper modules
import utils
import Spline

#import basic functionality
import Geometry
import DOS

#import modules which may not be present, depending on the distribution, ignore import errors
# path search modules
try:
	from Path import *
	import Calculators
	import Schedulers
except ImportError:
	pass

#  dimer method module
try:
	import Dimer
except ImportError:
	pass

__all__=["Geometry","Spline","constants","utils","DOS"]
