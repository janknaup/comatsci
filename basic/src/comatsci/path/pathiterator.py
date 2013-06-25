##############################################################################
# pathiterator.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2013 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


class PathIterator:
    """Iterate over Reactionpath instances"""
    
    def __init__(self,path):
        self.path=path
        self.index=0
    
    def __iter__(self):
        return self
    
    def next(self):
        if self.index==self.path.numimages():
            raise StopIteration
        else:
            self.index=self.index+1
            return self.path.geos[self.index]
