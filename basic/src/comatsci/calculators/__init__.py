##############################################################################
# Calculators/__init__.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

# define the list of known calculators here, also import base classes Calculator and CalcError
__all__=["noodlecalc", "siestacalc", "gaussiancalc", "erepcalc", 
          "muellerbrowncalc",  "calculator",  "calcerror",
          "pairpotentialcalc",  "potentials"]

# to maintain api compatibility with previous versions, import all calculator classes into the Calculators namespace
from noodlecalc  import *
from siestacalc  import *
from gaussiancalc  import *
from erepcalc   import *
from muellerbrowncalc  import *
from calculator import *
from calcerror import *
from pairpotentialcalc import pairPotentialCalc
import potentials

