##############################################################################
# ChargeConstraintsWizard.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################
# -*- coding: utf-8 -*-

from qt import *
from chargeConstraintsWizardBase import atomChargeConstraintsWizardBase


class atomChargeConstraintsWizard(atomChargeConstraintsWizardBase):

	def __init__(self,parent = None,name = None,modal = 0,fl = 0):
		atomChargeConstraintsWizardBase.__init__(self,parent,name,modal,fl)
		self.constraintsPrefactorLineEdit.setText("1.00")

	def pageFunctionDispatch(self, pageName):
		# call functions to be executed at page change
		if pageName=="Charge Constraints":
			self.printConstraints()


