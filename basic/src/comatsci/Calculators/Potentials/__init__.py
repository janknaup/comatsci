##############################################################################
# Calculators/Potentials/__init__.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

# define the list of known potentials here
__all__=["polynomFunction", "slaterFunction", "xySplineFunction", "potentialFunction" ]

# manually import potential functions

from potentialFunction import potentialFunction
from polynomFunction import polynomFunction
from slaterFunction import slaterFunction
from xySplineFunction import xySplineFunction


# define map of known potential functions
POTENTIALLABELMAP={
  "&polynom&": polynomFunction ,
  "&slater&": slaterFunction, 
  "&xyspline&"  :    xySplineFunction
}



def getPotentialFromFile(filename):
        """ read potential type from file and return respective potentialFunctionObject 
        @parameter filename: string containing input file name
        @return potentialFunctionObject from POTENTIALLABELMAP """
        #global POTENTIALLABELMAP
        infile=open(filename,"r")
        inlines=list(infile)
        infile.close()
        #first line should contain &POTENTIALLABEL& in the first field
        name=inlines[0].split()[0].strip()
        # check if first line actually contains a potential label
        if name[0]!="&" or name[-1]!="&":
            raise ValueError("Malformed potential file: Pair potential type indicator not found in line 1")
        else:
            # check if potentialFunction is known
            if not POTENTIALLABELMAP.has_key(name):
                raise ValueError("Unknown potential function '%s' specified" % name[1:-1])
            returnPotential= POTENTIALLABELMAP[name]({})
            returnPotential.parseString("\n".join(inlines))
        return returnPotential

