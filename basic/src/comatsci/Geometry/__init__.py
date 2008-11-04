##############################################################################
# Geometry/__init__.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

# import the basic functionality
from Geometry import *
from AnalysisGeometry import  *

#try to import Embedding support, ignore if not present in the current distribution
try:
	from EmbedGeometry import *
except ImportError:
	FullFeaturedGeometry=AnalysisGeometry
else:
	FullFeaturedGeometry=qmmmGeometry
