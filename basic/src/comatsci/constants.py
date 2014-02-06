##############################################################################
# constants.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

#
# This file contains various constants for global use in pastafarian/comatsci
#

#
# comatsci Version
#

from ._version import *  # @UnusedWildImport

# Verbosity levels:
# only fatal errors
VBL_SILENCE=-1
# fatal errors and warnings
VBL_QUIET=0
# VBL_QUIET+status and progress reports
VBL_NORMAL=1
#VBL_NORMAL+extended status reports
VBL_TALKY=2
#VBL_TALKY+basic debug information
VBL_DEBUG1=3
#VBL_DEBUG1+extended debug information
VBL_DEBUG2=4

#
# comatsci Calculators usable for pastafarian
#

PASTACALCS=['siesta','noodle','gaussian','muellerbrown']

#
# Unit conversion to atomic units
#
# This converts eV to H
EVOLT=0.036749325
#HARTREE=27.211384523 H to eV conversion
HARTREE=1.0

BOHR=1.0
ANGSTROM=5.291772e-01
NANOMETER=ANGSTROM/10.
METER=ANGSTROM*1e-10

KILOGRAM=1.660e-27

# This converts eV/Ang to H/Bohr
EVPERANG=0.01944691813482

#
# Miscellaneous stuff
#
#pi
PI=3.14159265358979
