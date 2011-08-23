# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'periodicExpand.ui'
#
# Created: Wed Aug 17 14:31:16 2011
#      by: The PyQt User Interface Compiler (pyuic) 3.18.1
#
# WARNING! All changes made in this file will be lost!


from qt import *
##############################################################################
# periodicExpand.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup , Knaup@bccms.uni-bremen.de
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################
#periodic expansion dialog


class periodicExpand(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QDialog.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("periodicExpand")

        self.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.sizePolicy().hasHeightForWidth()))
        self.setMinimumSize(QSize(300,310))
        self.setMaximumSize(QSize(300,310))
        self.setBaseSize(QSize(300,310))


        self.cDirGroupBox = QGroupBox(self,"cDirGroupBox")
        self.cDirGroupBox.setGeometry(QRect(10,190,280,60))

        self.cVectorLabel = QLabel(self.cDirGroupBox,"cVectorLabel")
        self.cVectorLabel.setGeometry(QRect(10,30,180,20))

        self.cSpinBox = QSpinBox(self.cDirGroupBox,"cSpinBox")
        self.cSpinBox.setGeometry(QRect(200,20,70,30))
        self.cSpinBox.setMinValue(1)

        self.textLabel1 = QLabel(self,"textLabel1")
        self.textLabel1.setGeometry(QRect(10,10,280,50))
        self.textLabel1.setAlignment(QLabel.WordBreak | QLabel.AlignVCenter)

        self.aDirGroupBox = QGroupBox(self,"aDirGroupBox")
        self.aDirGroupBox.setGeometry(QRect(10,70,280,60))

        self.aSpinBox = QSpinBox(self.aDirGroupBox,"aSpinBox")
        self.aSpinBox.setGeometry(QRect(200,20,70,30))
        self.aSpinBox.setMinValue(1)

        self.aVectorLabel = QLabel(self.aDirGroupBox,"aVectorLabel")
        self.aVectorLabel.setGeometry(QRect(10,30,180,20))

        self.cancelButton = QPushButton(self,"cancelButton")
        self.cancelButton.setGeometry(QRect(160,260,121,41))
        self.cancelButton.setAutoDefault(0)
        self.cancelButton.setDefault(1)

        self.bDirGroupBox = QGroupBox(self,"bDirGroupBox")
        self.bDirGroupBox.setGeometry(QRect(10,130,280,60))

        self.bVectorLabel = QLabel(self.bDirGroupBox,"bVectorLabel")
        self.bVectorLabel.setGeometry(QRect(10,30,180,20))

        self.bSpinBox = QSpinBox(self.bDirGroupBox,"bSpinBox")
        self.bSpinBox.setGeometry(QRect(200,20,70,30))
        self.bSpinBox.setMinValue(1)

        self.okButton = QPushButton(self,"okButton")
        self.okButton.setGeometry(QRect(20,260,120,41))
        self.okButton.setAutoDefault(0)
        self.okButton.setDefault(0)

        self.languageChange()

        self.resize(QSize(300,310).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.okButton,SIGNAL("clicked()"),self.accept)
        self.connect(self.cancelButton,SIGNAL("clicked()"),self.reject)

        self.setTabOrder(self.aSpinBox,self.bSpinBox)
        self.setTabOrder(self.bSpinBox,self.cSpinBox)
        self.setTabOrder(self.cSpinBox,self.cancelButton)
        self.setTabOrder(self.cancelButton,self.okButton)


    def languageChange(self):
        self.setCaption(self.__tr("Periodically Expand Geometry"))
        QToolTip.add(self,self.__tr("Geometry periodic expansion dialog"))
        QWhatsThis.add(self,self.__tr("This Dialog allows to periodically expand supercell geometries"))
        self.cDirGroupBox.setTitle(self.__tr("c-direction"))
        QToolTip.add(self.cDirGroupBox,self.__tr("c-vector"))
        QWhatsThis.add(self.cDirGroupBox,self.__tr("The x,y,z coordinates of the the supercell c-vector in Angstrom. Slect how many old supercells should be included in the new cell's c direction."))
        self.cVectorLabel.setText(QString.null)
        QToolTip.add(self.cSpinBox,self.__tr("number of old cells in c direction"))
        self.textLabel1.setText(self.__tr("<h2>Set number of old supercells in each lattice direction</h2>"))
        self.aDirGroupBox.setTitle(self.__tr("a-direction"))
        QToolTip.add(self.aDirGroupBox,self.__tr("a-vector"))
        QWhatsThis.add(self.aDirGroupBox,self.__tr("The x,y,z coordinates of the the supercell a-vector in Angstrom. Slect how many old supercells should be included in the new cell's a direction."))
        QToolTip.add(self.aSpinBox,self.__tr("number of old cells in a direction"))
        self.aVectorLabel.setText(QString.null)
        self.cancelButton.setText(self.__tr("cancel"))
        QToolTip.add(self.cancelButton,self.__tr("cancel periodic expansion"))
        self.bDirGroupBox.setTitle(self.__tr("b-direction"))
        QToolTip.add(self.bDirGroupBox,self.__tr("b-vector"))
        QWhatsThis.add(self.bDirGroupBox,self.__tr("The x,y,z coordinates of the the supercell b-vector in Angstrom. Slect how many old supercells should be included in the new cell's b direction."))
        self.bVectorLabel.setText(QString.null)
        QToolTip.add(self.bSpinBox,self.__tr("number of old cells in c direction"))
        self.okButton.setText(self.__tr("OK"))
        QToolTip.add(self.okButton,self.__tr("expand supercell and leave dialog"))


    def setVectors(self,a0):
        	Lattice=a0*5.291772e-01
        	self.aVectorLabel.setText("&lt;%7.4f %7.4f %7.4f&gt;" % (Lattice[0][0],Lattice[0][1],Lattice[0][2]))
        	self.bVectorLabel.setText("&lt;%7.4f %7.4f %7.4f&gt;" % (Lattice[1][0],Lattice[1][1],Lattice[1][2]))
        	self.cVectorLabel.setText("&lt;%7.4f %7.4f %7.4f&gt;" % (Lattice[2][0],Lattice[2][1],Lattice[2][2]))
        

    def getExpandSpecifier(self):
        	return (self.aSpinBox.value(),self.bSpinBox.value(),self.cSpinBox.value())
        

    def __tr(self,s,c = None):
        return qApp.translate("periodicExpand",s,c)
