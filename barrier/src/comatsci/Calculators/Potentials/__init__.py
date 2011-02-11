##############################################################################
# Calculators/Potentials/__init__.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

# define the list of known potentials here
__all__=["polynomFunction", "potentialFunction" ]

# manually import potential functions

from potentialFunction import potentialFunction
from polynomFunction import polynomFunction


# define map of known potential functions
POTENTIALLABELMAP={
  "&polynom&": polynomFunction 
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
          # try to construct potentialfunction from label, catch key errors from invalid potential labels and convert to ValueErrors
          try:
            returnPotential= POTENTIALLABELMAP[name]({})
          except KeyError as detail:
            raise ValueError("Unknown pair potential function %s" % detail.args[0][1:-1])
        returnPotential.parseString("\n".join(inlines))
        return returnPotential

