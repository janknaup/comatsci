##############################################################################
# <filename>
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################
# -*- coding: utf-8 -*-

from qt import *
from linkLayersBase import layerLinkBase


class layerLinkHCS(layerLinkBase):

	def __init__(self,parent = None,name = None,modal = 0,fl = 0):
		layerLinkBase.__init__(self,parent,name,modal,fl)
	
	# public slot
	def linkAction(self):
		# get the link-atom distance factor
		ladf=float(self.LADFlineEdit.text())
		# create a progress dialog and execute the linking operation
		progressWindow=QProgressDialog()
		progressWindow.setLabelText("generating HCS link atoms")
		progressWindow.show()
		(self.embeddedGeometry,results)=self.geo.layersubgeometry(self.QMZComboBox.currentItem()).HCSLinkedGeometry(self.geo.layersubgeometry(self.PCHRComboBox.currentItem()),progressWindow.setTotalSteps,progressWindow.setProgress,distscale=ladf)
		progressWindow.hide()
		return results
