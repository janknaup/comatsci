# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'chargeConstraintsWizardBase.ui'
#
# Created: Wed Aug 17 14:31:16 2011
#      by: The PyQt User Interface Compiler (pyuic) 3.18.1
#
# WARNING! All changes made in this file will be lost!


from qt import *
##############################################################################
# chargeConstraintsWizardBase.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup, Knaup@bccms.uni-bremen.de
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


class atomChargeConstraintsWizardBase(QWizard):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QWizard.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("chargeConstraintsWizard")



        self.page = QWidget(self,"page")
        pageLayout = QGridLayout(self.page,1,1,6,6,"pageLayout")

        self.textLabel1 = QLabel(self.page,"textLabel1")
        self.textLabel1.setTextFormat(QLabel.AutoText)
        self.textLabel1.setAlignment(QLabel.WordBreak | QLabel.AlignVCenter)

        pageLayout.addMultiCellWidget(self.textLabel1,0,0,0,1)

        self.atomSerialsTextEdit = QTextEdit(self.page,"atomSerialsTextEdit")
        self.atomSerialsTextEdit.setTextFormat(QTextEdit.PlainText)

        pageLayout.addMultiCellWidget(self.atomSerialsTextEdit,2,2,0,1)

        self.splitter1 = QSplitter(self.page,"splitter1")
        self.splitter1.setOrientation(QSplitter.Horizontal)

        self.textLabel1_2 = QLabel(self.splitter1,"textLabel1_2")

        self.constraintsPrefactorLineEdit = QLineEdit(self.splitter1,"constraintsPrefactorLineEdit")

        pageLayout.addWidget(self.splitter1,1,0)
        spacer1 = QSpacerItem(151,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
        pageLayout.addItem(spacer1,1,1)
        self.addPage(self.page,QString(""))

        self.WizardPage = QWidget(self,"WizardPage")
        WizardPageLayout = QGridLayout(self.WizardPage,1,1,6,6,"WizardPageLayout")

        self.constraintsDisplayTextEdit = QTextEdit(self.WizardPage,"constraintsDisplayTextEdit")
        self.constraintsDisplayTextEdit.setReadOnly(1)

        WizardPageLayout.addWidget(self.constraintsDisplayTextEdit,0,0)
        self.addPage(self.WizardPage,QString(""))

        self.languageChange()

        self.resize(QSize(600,480).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self,SIGNAL("selected(const QString&)"),self.pageFunctionDispatch)


    def languageChange(self):
        self.setCaption(self.__tr("DFTB+ Charge Contraints"))
        QToolTip.add(self,self.__tr("Generate atom charge constraint specifications in DFTB+ .hsd format."))
        self.textLabel1.setText(self.__tr("<font size=\"+1\">Input serial numbers (counting from 1) of atoms to apply constraints on, separated by white space, in the text editor below. <br>When ready, click \"next\".</font>"))
        QToolTip.add(self.atomSerialsTextEdit,self.__tr("Enter the serial numbers of atoms here. Atom serials can i.e. be obtained from vmd by using atomselection get serial."))
        self.textLabel1_2.setText(self.__tr("charge constraints prefactor"))
        self.constraintsPrefactorLineEdit.setText(self.__tr("."))
        self.constraintsPrefactorLineEdit.setInputMask(self.__tr("D.dd;0"))
        self.setTitle(self.page,self.__tr("Select Atoms"))
        QToolTip.add(self.constraintsDisplayTextEdit,self.__tr("Copy and paste constraints specifiers into DFTB+ dftb_in.hsd file. See DFTB+ manual!"))
        self.setTitle(self.WizardPage,self.__tr("Charge Constraints"))


    def printConstraints(self):
        	# parse the atom serial numbers, get the DFTB+ charge constraints string, put them into the output textEdit and switch pages
        	serials=str(self.atomSerialsTextEdit.text())
        	# user interface counts from 1, internal indices count from 0
        	try:
        		indices=[ int(j)-1 for j in serials.split() ]
        	except ValueError:
        		QMessageBox.warning(self,"Invalid atom numbers","<H1>Error</H1><p>The list of atom serials could not be understood. Atom serials must be given as a white-space separated list of integers. Please check your input.</p>")
        	# convert prefactor to float
        	prefactor=float(self.constraintsPrefactorLineEdit.text())
        	try:
        		self.constraintsDisplayTextEdit.setText(self.geo.getHSDChargeConstraintString(indices, prefactor))
        	except ValueError:
        		QMessageBox.warning(self,"Invalid atom number","<H1>Error</H1><p>One or more atom serials are not present in the geometry. Please check your input.</p>")
        	except:
        		raise
        	# enable finish button
        	self.setFinishEnabled(self.WizardPage,True)
        

    def pageFunctionDispatch(self):
        print "atomChargeConstraintsWizardBase.pageFunctionDispatch(): Not implemented yet"

    def setGeometry(self,a0):
        	# store a link to the geometry object to work on
        	self.geo=a0
        

    def __tr(self,s,c = None):
        return qApp.translate("atomChargeConstraintsWizardBase",s,c)
