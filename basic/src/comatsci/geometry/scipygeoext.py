##############################################################################
# scipygeoext.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2013 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################
#
# drop-in replacement for old c-extensions using scipy 0.12 or higher
#
##############################################################################

import numpy
from scipy import spatial
from .. import constants

def dmatrix(geo):
    """return the (n,n) distance matrix for a numarray of (n,3) float coordinates
    @param geo: ndarray
    @param geo: geometry array
    @return: (n,n) ndarray of all distances between points in geo"""
    # almost too trivial
    return spatial.distance.cdist(geo,geo,metric='euclidean')


def sdmatrix(geo,lattice):
    """return the (n,n) minimum distance matrix for a numarray of (n,3) float coordinates with (3,3) array of lattice vectors
    @type geo: ndarray
    @param geo: geometry array
    @type lattice: (3,3) ndarray
    @param lattice: periodic boundary translation vectors
    @return: (n,n) ndarray of all minimum distances between points in geo and periodic images
    @see: Calls the L{crossSupercellDistanceMatrix} function"""
    #see below
    return crossSupercellDistanceMatrix(geo,geo,lattice)


def crossSupercellDistanceMatrix(a,b,lattice):
    """return the (n,m) minimum distance matrix for two numarrays of (n,3) and (m,3) float coordinates with (3,3) array of lattice vectors
    @type a: ndarray
    @param a: geometry array
    @type b: ndarray
    @param b: geometry array
    @type lattice: (3,3) ndarray
    @param lattice: periodic boundary translation vectors
    @return: (n,n) ndarray of all minimum distances between points in geo and periodic images
    @warning: uses (2*)27*n*n floats temporary memory, this may be a problem with large geometries
    """
    # this could be rewritten to use about 1/10 for the price of more CPU calls by comparing
    # distances to each periodic image directly after calculating them instead of accumulating first 
    dms=[]  # distance matrices to the periodic images are stored here
    # iterate over one neighbor image per direction.
    for i in range(-1,2,1):
        for j in range(-1,2,1):
            for k in range (-1,2,1):
                dms.append(spatial.distance.cdist(a+(lattice*numpy.array((i,j,k))).sum(axis=1),b,metric='euclidean'))
    return numpy.array(dms).min(axis=0)
    
    
def blmatrix(types,corad):
    """return the (n,n) matrix of covalent radius sums for the types and covalent radius lists provided
    @type types: iterable of integers
    @param types: atom types list (Z-numbers)
    @type corad: sequence of floats
    @param corad: sequence of covalent radii in order of Z
    @rtype: ndarray of floats
    @return: matrix of covalent bond distances between atom types listed in types     
    """
    ARAD=numpy.array(corad)[types]
    return numpy.array(numpy.meshgrid(ARAD,ARAD)).sum(axis=0)


def blist(positions,types,corad,tolerance=1.1):
    """return the per-atom list of bond parter lists, molecule version
    @warning: Expects positions in Bohr and corad in Angstrom! Stupid weird behavior kept to stay drop-in
        compatible. May be corrected in the future.
    @type positions: (n,3) ndarray of floats
    @param positions: the tom positions to find bonds between
    @type types: sequence of integer
    @param types: types list of the atoms to search bonds between 
    @type corad: sequence of floats
    @param corad: sequence of covalent radii in order of Z
    @rtype: list of lists of integers
    @return: list of lists containing bond partners per atom in sequence
    @see: Calls the L{blmatrix} and L{dmatrix} functions"""
    dm=dmatrix(positions)
    BB=(blmatrix(types,corad)/constants.ANGSTROM)*tolerance
    bondmatrix=numpy.less(dm,BB*(numpy.ones_like(BB)-numpy.eye(BB.shape[0]))*1.1)
    bondlist=[]
    for ii in bondmatrix:
        bondlist.append(numpy.nonzero(ii)[0].tolist())
    return bondlist
    


def sblist(positions,lattice,types,corad,tolerance=1.1):
    """return a tuple of the per-atom list of bond parter lists and bond partner periodic image coordinates, supercell version
    @type positions: (n,3) ndarray of floats
    @param positions: the tom positions to find bonds between
    @type lattice: (3,3) ndarray
    @param lattice: periodic boundary translation vectors
    @type types: sequence of integer
    @param types: types list of the atoms to search bonds between 
    @type corad: sequence of floats
    @param corad: sequence of covalent radii in order of Z
    @rtype: list of lists of integers
    @return: list of lists containing bond partners per atom in sequence
    @see: Calls the L{blmatrix} and L{dmatrix} functions. sdmatrix and crossSupercellDistanceMatrix are B{not} called."""
    BB=(blmatrix(types,corad)/constants.ANGSTROM)*tolerance
    natom=len(types)
    bondlist=[[]]*natom
    imagecoordlist=[[]]*natom
    for a in range(-1,2,1):
        for b in range(-1,2,1):
            for c in range (-1,2,1):
                sccoord=[(a,b,c)]
                dm=spatial.distance.cdist(positions+(lattice*numpy.array((a,b,c))).sum(axis=1),positions)
                bondmatrix=numpy.less(dm,BB*(numpy.ones_like(BB)-numpy.eye(BB.shape[0]))*1.1)
                for ii in range(natom):
                    tempbonds=numpy.nonzero(bondmatrix[ii])[0].tolist()
                    bondlist[ii]=bondlist[ii]+tempbonds
                    imagecoordlist[ii]+=sccoord*len(tempbonds)
    for ii in range(natom):
        bondlist[ii].sort()
    return [bondlist,imagecoordlist]
