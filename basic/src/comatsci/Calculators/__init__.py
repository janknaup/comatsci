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
__all__=["dftbcalc", "noodlecalc", "siestacalc", "gaussiancalc", "erepcalc",  "muellerBrownCalc",  "Calculator",  "CalcError",
          "pairPotentialCalc",  "Potentials"]

# to maintain api compatibility with previous versions, import all calculatir classes into the Calculators namespace
# @note: old-version dftb calculator removed from official distribution. noodle is now the officially supported DFTB version
# from dftbcalc  import *
from noodlecalc  import *
from siestacalc  import *
from gaussiancalc  import *
from erepcalc   import *
from muellerBrownCalc  import *
from Calculator   import *
from CalcError import *
from pairPotentialCalc import pairPotentialCalc
import Potentials

