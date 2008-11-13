# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'EditGeometry.ui'
#
# Created: Thu Nov 13 11:18:05 2008
#      by: The PyQt User Interface Compiler (pyuic) 3.17.4
#
# WARNING! All changes made in this file will be lost!


from qt import *
from qttable import QTable
##############################################################################
# EditGeometry.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup , Knaup@bccms.uni-bremen.de
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

#Atom editing form
from comatsci import Geometry
from qttable import QTableItem
from qttable import QComboTableItem
import copy


class Edit_Atoms(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QDialog.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("Edit_Atoms")

        self.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.sizePolicy().hasHeightForWidth()))
        self.setMinimumSize(QSize(780,430))
        self.setMaximumSize(QSize(780,430))
        self.setBaseSize(QSize(780,430))


        self.frame6 = QFrame(self,"frame6")
        self.frame6.setGeometry(QRect(530,340,250,90))
        self.frame6.setFrameShape(QFrame.StyledPanel)
        self.frame6.setFrameShadow(QFrame.Raised)

        self.okButton = QPushButton(self.frame6,"okButton")
        self.okButton.setGeometry(QRect(10,10,110,31))

        self.revertButton = QPushButton(self.frame6,"revertButton")
        self.revertButton.setGeometry(QRect(130,10,110,31))

        self.cancelButton = QPushButton(self.frame6,"cancelButton")
        self.cancelButton.setGeometry(QRect(130,50,110,31))

        self.applyButton = QPushButton(self.frame6,"applyButton")
        self.applyButton.setGeometry(QRect(10,50,110,31))

        self.frame5 = QFrame(self,"frame5")
        self.frame5.setGeometry(QRect(0,280,520,150))
        self.frame5.setFrameShape(QFrame.StyledPanel)
        self.frame5.setFrameShadow(QFrame.Raised)

        self.layerNameEdit = QLineEdit(self.frame5,"layerNameEdit")
        self.layerNameEdit.setGeometry(QRect(10,60,210,31))

        self.layersFrame = QLabel(self.frame5,"layersFrame")
        self.layersFrame.setGeometry(QRect(10,10,101,33))

        self.layerTable = QTable(self.frame5,"layerTable")
        self.layerTable.setNumCols(self.layerTable.numCols() + 1)
        self.layerTable.horizontalHeader().setLabel(self.layerTable.numCols() - 1,self.__tr("#"))
        self.layerTable.setNumCols(self.layerTable.numCols() + 1)
        self.layerTable.horizontalHeader().setLabel(self.layerTable.numCols() - 1,self.__tr("Name"))
        self.layerTable.setGeometry(QRect(230,10,280,130))
        self.layerTable.setNumRows(3)
        self.layerTable.setNumCols(2)

        self.newLayerButton = QPushButton(self.frame5,"newLayerButton")
        self.newLayerButton.setGeometry(QRect(9,100,100,31))

        self.deleteLayerButton = QPushButton(self.frame5,"deleteLayerButton")
        self.deleteLayerButton.setGeometry(QRect(120,100,100,31))

        self.frame7 = QFrame(self,"frame7")
        self.frame7.setGeometry(QRect(530,280,250,50))
        self.frame7.setFrameShape(QFrame.StyledPanel)
        self.frame7.setFrameShadow(QFrame.Raised)

        self.addAtomButton = QPushButton(self.frame7,"addAtomButton")
        self.addAtomButton.setGeometry(QRect(10,10,111,31))

        self.deleteAtomButton = QPushButton(self.frame7,"deleteAtomButton")
        self.deleteAtomButton.setGeometry(QRect(130,10,110,31))

        self.atomTable = QTable(self,"atomTable")
        self.atomTable.setNumCols(self.atomTable.numCols() + 1)
        self.atomTable.horizontalHeader().setLabel(self.atomTable.numCols() - 1,self.__tr("elem"))
        self.atomTable.setNumCols(self.atomTable.numCols() + 1)
        self.atomTable.horizontalHeader().setLabel(self.atomTable.numCols() - 1,self.__tr("x"))
        self.atomTable.setNumCols(self.atomTable.numCols() + 1)
        self.atomTable.horizontalHeader().setLabel(self.atomTable.numCols() - 1,self.__tr("y"))
        self.atomTable.setNumCols(self.atomTable.numCols() + 1)
        self.atomTable.horizontalHeader().setLabel(self.atomTable.numCols() - 1,self.__tr("z"))
        self.atomTable.setNumCols(self.atomTable.numCols() + 1)
        self.atomTable.horizontalHeader().setLabel(self.atomTable.numCols() - 1,self.__tr("layer"))
        self.atomTable.setNumCols(self.atomTable.numCols() + 1)
        self.atomTable.horizontalHeader().setLabel(self.atomTable.numCols() - 1,self.__tr("charge"))
        self.atomTable.setNumCols(self.atomTable.numCols() + 1)
        self.atomTable.horizontalHeader().setLabel(self.atomTable.numCols() - 1,self.__tr("subtype"))
        self.atomTable.setGeometry(QRect(0,0,780,270))
        self.atomTable.setNumRows(12)
        self.atomTable.setNumCols(7)
        self.atomTable.setReadOnly(0)
        self.atomTable.setSelectionMode(QTable.SingleRow)

        self.languageChange()

        self.resize(QSize(780,430).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.newLayerButton,SIGNAL("clicked()"),self.addLayer)
        self.connect(self.revertButton,SIGNAL("clicked()"),self.revert)
        self.connect(self.cancelButton,SIGNAL("clicked()"),self.close)
        self.connect(self.applyButton,SIGNAL("clicked()"),self.my_apply)
        self.connect(self.okButton,SIGNAL("clicked()"),self.my_ok)
        self.connect(self.deleteLayerButton,SIGNAL("clicked()"),self.deleteLayer)
        self.connect(self.addAtomButton,SIGNAL("clicked()"),self.addAtom)
        self.connect(self.deleteAtomButton,SIGNAL("clicked()"),self.deleteAtom)

        self.setTabOrder(self.cancelButton,self.revertButton)
        self.setTabOrder(self.revertButton,self.applyButton)
        self.setTabOrder(self.applyButton,self.okButton)
        self.setTabOrder(self.okButton,self.addAtomButton)
        self.setTabOrder(self.addAtomButton,self.deleteAtomButton)
        self.setTabOrder(self.deleteAtomButton,self.newLayerButton)
        self.setTabOrder(self.newLayerButton,self.deleteLayerButton)
        self.setTabOrder(self.deleteLayerButton,self.layerNameEdit)
        self.setTabOrder(self.layerNameEdit,self.layerTable)
        self.setTabOrder(self.layerTable,self.atomTable)


    def languageChange(self):
        self.setCaption(self.__tr("Edit Atoms"))
        QToolTip.add(self,self.__tr("Atom Editor Dialog"))
        QWhatsThis.add(self,self.__tr("This dialog allows to edit the atoms and layers of a geometry"))
        QWhatsThis.add(self.frame6,self.__tr("Accept, apply or revert changes or cancel"))
        self.okButton.setText(self.__tr("OK"))
        QToolTip.add(self.okButton,self.__tr("accept changes and exit dialog"))
        self.revertButton.setText(self.__tr("revert"))
        QToolTip.add(self.revertButton,self.__tr("revert all changes since last apply"))
        self.cancelButton.setText(self.__tr("cancel"))
        QToolTip.add(self.cancelButton,self.__tr("revert changes and exit dialog"))
        self.applyButton.setText(self.__tr("apply"))
        QToolTip.add(self.applyButton,self.__tr("apply changes"))
        QToolTip.add(self.frame5,self.__tr("layer editor"))
        QWhatsThis.add(self.frame5,self.__tr("add or delete layers to the geometry"))
        self.layersFrame.setText(self.__tr("<h1>Layers</h1>"))
        self.layerTable.horizontalHeader().setLabel(0,self.__tr("#"))
        self.layerTable.horizontalHeader().setLabel(1,self.__tr("Name"))
        QToolTip.add(self.layerTable,self.__tr("Table of layers"))
        QWhatsThis.add(self.layerTable,self.__tr("This Table lists all layers in the geometry"))
        self.newLayerButton.setText(self.__tr("add new"))
        QToolTip.add(self.newLayerButton,self.__tr("append a new layer to the list"))
        self.deleteLayerButton.setText(self.__tr("delete"))
        QToolTip.add(self.deleteLayerButton,self.__tr("delete the currently selected layer"))
        QWhatsThis.add(self.frame7,self.__tr("add or delete atoms from the current geometry"))
        self.addAtomButton.setText(self.__tr("add atom"))
        QToolTip.add(self.addAtomButton,self.__tr("add a new atom before the currently selected"))
        self.deleteAtomButton.setText(self.__tr("delete atom"))
        QToolTip.add(self.deleteAtomButton,self.__tr("delete the currently selected atom"))
        self.atomTable.horizontalHeader().setLabel(0,self.__tr("elem"))
        self.atomTable.horizontalHeader().setLabel(1,self.__tr("x"))
        self.atomTable.horizontalHeader().setLabel(2,self.__tr("y"))
        self.atomTable.horizontalHeader().setLabel(3,self.__tr("z"))
        self.atomTable.horizontalHeader().setLabel(4,self.__tr("layer"))
        self.atomTable.horizontalHeader().setLabel(5,self.__tr("charge"))
        self.atomTable.horizontalHeader().setLabel(6,self.__tr("subtype"))
        QToolTip.add(self.atomTable,self.__tr("Table of atoms"))
        QWhatsThis.add(self.atomTable,self.__tr("This table lists all atoms in the geometry"))


    def addLayer(self):
        	lname=str(self.layerNameEdit.text())
        	li=self.geo.addlayer(lname)
        	newrow=self.layerTable.numRows()
        	self.layerTable.setNumRows(newrow+1)
        	self.layerTable.setItem(newrow,0,QTableItem(self.layerTable,QTableItem.Never,str(li)))
        	self.layerTable.setItem(newrow,1,QTableItem(self.layerTable,QTableItem.Never,lname))
        	layerlabelstring="%i: %s" % (li,lname)
        	self.layersqsl.append(layerlabelstring)
        	self.updateLayersInAtomTable()
        	self.layerNameEdit.setText("")
        		

    def my_ok(self):
        	self.my_apply()
        	self.close()
        

    def setGeometry(self,a0):
        	#keep a the original geometry untouched
        	self.originalgeo=a0
        	self.geo=copy.deepcopy(a0)
        	#initialize layer table and Layer Spinboxes
        	layercount=len(self.geo.LayerDict)
        	self.layerTable.setNumRows(layercount)
        	self.layersqsl=QStringList()
        	for i in range(layercount):
        		li=self.geo.LayerDict.keys()[i]
        		lname=self.geo.LayerDict[li].Name
        		layerlabelstring="%i: %s" % (li,lname)
        		self.layerTable.setItem(i,0,QTableItem(self.layerTable,QTableItem.Never,str(li)))
        		self.layerTable.setItem(i,1,QTableItem(self.layerTable,QTableItem.Never,str(lname)))
        		self.layersqsl.append(layerlabelstring)
        	self.layerTable.adjustColumn(0)
        	self.layerTable.setColumnStretchable(1,True)
        	#initialize atom table
        	self.PTEqsl=QStringList()
        	for i in self.geo.PTE:
        		self.PTEqsl.append(i)
        	self.atomTable.setNumRows(self.geo.Atomcount)
        	for i in range(self.geo.Atomcount):
        		self.atomTable.setItem(i,0,QComboTableItem(self.atomTable,self.PTEqsl,False))
        		self.atomTable.item(i,0).setCurrentItem(self.geo.AtomTypes[i])
        		self.atomTable.setItem(i,1,QTableItem(self.atomTable,QTableItem.WhenCurrent,str(self.geo.Geometry[i][0])))
        		self.atomTable.setItem(i,2,QTableItem(self.atomTable,QTableItem.WhenCurrent,str(self.geo.Geometry[i][1])))
        		self.atomTable.setItem(i,3,QTableItem(self.atomTable,QTableItem.WhenCurrent,str(self.geo.Geometry[i][2])))
        		self.atomTable.setItem(i,4,QComboTableItem(self.atomTable,self.layersqsl,False))
        		self.atomTable.item(i,4).setCurrentItem(self.geo.AtomLayers[i])
        		self.atomTable.setItem(i,5,QTableItem(self.atomTable,QTableItem.WhenCurrent,str(self.geo.AtomCharges[i])))
        		self.atomTable.setItem(i,6,QTableItem(self.atomTable,QTableItem.WhenCurrent,str(self.geo.AtomSubTypes[i])))
        	for i in range(6):
        		self.atomTable.adjustColumn(i)
        	self.atomTable.setColumnStretchable(6,True)
        	

    def deleteAtom(self):
        	curatom=self.atomTable.currentRow()
        	self.atomTable.removeRow(curatom)
        

    def revert(self):
        	self.setGeometry(self.originalgeo)
        

    def my_apply(self):
        	#construct a  completely new geometry object for output
        	newgeo=self.geo.__class__()
        	#first set the new geometry's Layer Dictionary     '
        	for i in range(self.layerTable.numRows()):
        		li=int(str(self.layerTable.text(i,0)))
        		lname=str(self.layerTable.text(i,1))
        		newgeo.LayerDict[li]=Geometry.GeoLayer(lname)
        	#just reuse mode Lattice and origin
        	newgeo.Origin=self.geo.Origin
        	newgeo.Lattice=self.geo.Lattice
        	newgeo.Mode=self.geo.Mode
        	#now populate the new Geometry with atoms
        	for i in range(self.atomTable.numRows()):
        		type=self.atomTable.item(i,0).currentItem()
        		x=float(str(self.atomTable.text(i,1)))
        		y=float(str(self.atomTable.text(i,2)))
        		z=float(str(self.atomTable.text(i,3)))
        		li=self.atomTable.item(i,4).currentItem()
        		chr=float(str(self.atomTable.text(i,5)))
        		st=str(self.atomTable.text(i,6))
        		newgeo.addatom(type,(x,y,z),li,chr,st)
        	#This cunningly overwrites the original geometry object
        	self.setGeometry(newgeo)
        

    def addAtom(self):
        	curatom=self.atomTable.currentRow()
        	type=self.atomTable.item(curatom,0).currentItem()
        	li=self.atomTable.item(curatom,4).currentItem()
        	chr=self.atomTable.text(curatom,5)
        	st=self.atomTable.text(curatom,6)
        	self.atomTable.insertRows(curatom)
        	self.atomTable.setItem(curatom,0,QComboTableItem(self.atomTable,self.PTEqsl,False))
        	self.atomTable.item(curatom,0).setCurrentItem(type)
        	self.atomTable.setText(curatom,1,"")
        	self.atomTable.setText(curatom,2,"")
        	self.atomTable.setText(curatom,3,"")
        	self.atomTable.setItem(curatom,4,QComboTableItem(self.atomTable,self.layersqsl,False))
        	self.atomTable.item(curatom,4).setCurrentItem(li)
        	self.atomTable.setText(curatom,5,"0.0")
        	self.atomTable.setText(curatom,6,st)
        

    def updateLayersInAtomTable(self):
        	for i in range(self.atomTable.numRows()):
        		olditem=self.atomTable.item(i,4).currentItem()
        		self.atomTable.setItem(i,4,QComboTableItem(self.atomTable,self.layersqsl,False))
        		self.atomTable.item(i,4).setCurrentItem(olditem)
        

    def deleteLayer(self):
        	curlayer=self.layerTable.currentRow()
        	# Do not touch the default layer
        	if curlayer!=0:
        		self.layerTable.removeRow(curlayer)
        		#renumber the layers
        		for i in range(curlayer-1,self.layerTable.numRows()):
        			self.layerTable.setText(i,0,str(i))
        		#put atoms from deleted layer into default layer
        		#and make atoms with layer indices higher than deleted layer consistent again
        		for i in range(self.atomTable.numRows()):
        			itmlr=self.atomTable.item(i,4).currentItem()
        			if itmlr==curlayer:
        				self.atomTable.item(i,4).setCurrentItem(0)
        			elif itmlr>curlayer:
        				self.atomTable.item(i,4).setCurrentItem(itmlr-1)
        		#make layers QStringList consistent with table and update comboBoxes for all atoms
        		del self.layersqsl[curlayer]
        		self.updateLayersInAtomTable()
        

    def __tr(self,s,c = None):
        return qApp.translate("Edit_Atoms",s,c)
