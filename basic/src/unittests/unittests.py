'''
Created on Mar 13, 2012

@author: Jan M. Knaup <janknaup@gmail.com>
'''

from __future__ import print_function

import unittest

testNames=("testGeometry",
          "testSpline")

tests=map(__import__, testNames)

class comatsciTest(unittest.TestSuite):
    
    def __init__(self):
        unittest.TestSuite.__init__(self,tests)


if __name__ == '__main__':
    tests=map(__import__, testNames)
    print(tests)
    allsuite=unittest.TestSuite(tests)
    unittest.TextTestRunner(verbosity=2).run(allsuite)