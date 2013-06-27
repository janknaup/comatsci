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
        if self.index==self.path.numimages()-1:
            raise StopIteration
        else:
            self.index=self.index+1
            return self.path.geos[self.index]


class EnergyAccessor:
    def __init__(self,path):
        self.path=path
        self.index=0
        
    def __iter__(self):
        return self
    
    def next(self):
        if self.index==self.path.numimages()-1:
            raise StopIteration
        else:
            self.index=self.index+1
            return self.path.geos[self.index].totalenergy
    
    def __getitem__(self,index):
        if isinstance(index,slice):
            if index.start==None: 
                start=0
            else:
                start=index.start
            if index.step==None:
                step=1
            else:
                step=index.step
            return [ self.__getitem__(ii) for ii in range(start,index.stop,step) ]
        else:
            if index<0: index=index+self.__len__()
            if index < 0 or index>=self.path.numImages:
                raise IndexError
            else:
                return self.path.geos[index].totalenergy
    
    def __setitem__(self,index,value):
        if abs(index)>=self.path.numImages:
            raise IndexError
        else:
            self.path.geos[index].totalenergy=value
            
    def __len__(self):
        return len(self.path.geos)


class ForcesAccessor:
    def __init__(self,path):
        self.path=path
        self.index=0
        
    def __iter__(self):
        return self
    
    def next(self):
        if self.index==self.path.numimages()-1:
            raise StopIteration
        else:
            self.index=self.index+1
            return self.path.geos[self.index].forces
    
    def __getitem__(self,index):
        if isinstance(index,slice):
            if index.start==None: 
                start=0
            else:
                start=index.start
            if index.step==None:
                step=1
            else:
                step=index.step
            return [ self.__getitem__(ii) for ii in range(start,index.stop,step) ]
        else:
            if index<0: index=index+self.__len__()
            if index < 0 or index>=self.path.numImages:
                raise IndexError
            else:
                return self.path.geos[index].forces
    
    def __setitem__(self,index,value):
        if abs(index)>=self.path.numImages:
            raise IndexError
        else:
            self.path.geos[index].forces=value
            
    def __len__(self):
        return len(self.path.geos)