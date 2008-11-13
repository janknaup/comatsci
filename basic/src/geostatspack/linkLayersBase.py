# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'linkLayersBase.ui'
#
# Created: Thu Nov 13 11:18:06 2008
#      by: The PyQt User Interface Compiler (pyuic) 3.17.4
#
# WARNING! All changes made in this file will be lost!


from qt import *
##############################################################################
# linkLayersBase.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup , Knaup@bccms.uni-bremen.de
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from comatsci import utils
import copy


class layerLinkBase(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QDialog.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("LayerLinkBase")

        self.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Minimum,0,0,self.sizePolicy().hasHeightForWidth()))
        self.setMinimumSize(QSize(600,400))
        self.setBaseSize(QSize(600,450))
        self.setModal(1)

        LayerLinkBaseLayout = QVBoxLayout(self,6,6,"LayerLinkBaseLayout")

        self.groupBox1 = QGroupBox(self,"groupBox1")
        self.groupBox1.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Preferred,0,0,self.groupBox1.sizePolicy().hasHeightForWidth()))
        self.groupBox1.setMinimumSize(QSize(590,150))
        self.groupBox1.setColumnLayout(0,Qt.Vertical)
        self.groupBox1.layout().setSpacing(6)
        self.groupBox1.layout().setMargin(6)
        groupBox1Layout = QVBoxLayout(self.groupBox1.layout())
        groupBox1Layout.setAlignment(Qt.AlignTop)

        self.frame5 = QFrame(self.groupBox1,"frame5")
        self.frame5.setFrameShape(QFrame.NoFrame)
        self.frame5.setFrameShadow(QFrame.Raised)
        frame5Layout = QHBoxLayout(self.frame5,6,6,"frame5Layout")

        self.textLabel1 = QLabel(self.frame5,"textLabel1")
        frame5Layout.addWidget(self.textLabel1)

        self.LADFlineEdit = QLineEdit(self.frame5,"LADFlineEdit")
        frame5Layout.addWidget(self.LADFlineEdit)
        spacer2 = QSpacerItem(350,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
        frame5Layout.addItem(spacer2)
        groupBox1Layout.addWidget(self.frame5)

        self.frame6 = QFrame(self.groupBox1,"frame6")
        self.frame6.setFrameShape(QFrame.NoFrame)
        self.frame6.setFrameShadow(QFrame.Raised)
        frame6Layout = QGridLayout(self.frame6,1,1,6,6,"frame6Layout")

        self.frame7 = QFrame(self.frame6,"frame7")
        self.frame7.setFrameShape(QFrame.NoFrame)
        self.frame7.setFrameShadow(QFrame.Raised)
        frame7Layout = QGridLayout(self.frame7,1,1,6,6,"frame7Layout")

        self.textLabel2 = QLabel(self.frame7,"textLabel2")

        frame7Layout.addWidget(self.textLabel2,0,1)

        self.QMZComboBox = QComboBox(0,self.frame7,"QMZComboBox")

        frame7Layout.addMultiCellWidget(self.QMZComboBox,1,1,0,2)

        frame6Layout.addWidget(self.frame7,0,0)

        self.frame7_2 = QFrame(self.frame6,"frame7_2")
        self.frame7_2.setFrameShape(QFrame.NoFrame)
        self.frame7_2.setFrameShadow(QFrame.Raised)
        frame7_2Layout = QGridLayout(self.frame7_2,1,1,6,6,"frame7_2Layout")

        self.textLabel2_2 = QLabel(self.frame7_2,"textLabel2_2")

        frame7_2Layout.addWidget(self.textLabel2_2,0,1)

        self.PCHRComboBox = QComboBox(0,self.frame7_2,"PCHRComboBox")

        frame7_2Layout.addMultiCellWidget(self.PCHRComboBox,1,1,0,2)

        frame6Layout.addWidget(self.frame7_2,0,1)
        groupBox1Layout.addWidget(self.frame6)
        LayerLinkBaseLayout.addWidget(self.groupBox1)

        self.groupBox2 = QGroupBox(self,"groupBox2")
        self.groupBox2.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Expanding,0,0,self.groupBox2.sizePolicy().hasHeightForWidth()))
        self.groupBox2.setMinimumSize(QSize(590,200))
        self.groupBox2.setColumnLayout(0,Qt.Vertical)
        self.groupBox2.layout().setSpacing(6)
        self.groupBox2.layout().setMargin(6)
        groupBox2Layout = QGridLayout(self.groupBox2.layout())
        groupBox2Layout.setAlignment(Qt.AlignTop)

        self.resultsTextEdit = QTextEdit(self.groupBox2,"resultsTextEdit")
        self.resultsTextEdit.setEnabled(1)
        self.resultsTextEdit.setTextFormat(QTextEdit.PlainText)
        self.resultsTextEdit.setWordWrap(QTextEdit.NoWrap)
        self.resultsTextEdit.setReadOnly(1)
        self.resultsTextEdit.setUndoRedoEnabled(0)

        groupBox2Layout.addWidget(self.resultsTextEdit,0,0)
        LayerLinkBaseLayout.addWidget(self.groupBox2)

        self.buttonsFrame = QFrame(self,"buttonsFrame")
        self.buttonsFrame.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.buttonsFrame.sizePolicy().hasHeightForWidth()))
        self.buttonsFrame.setMinimumSize(QSize(590,50))
        self.buttonsFrame.setFrameShape(QFrame.StyledPanel)
        self.buttonsFrame.setFrameShadow(QFrame.Raised)
        buttonsFrameLayout = QHBoxLayout(self.buttonsFrame,6,6,"buttonsFrameLayout")

        self.linkButton = QPushButton(self.buttonsFrame,"linkButton")
        buttonsFrameLayout.addWidget(self.linkButton)
        spacer1 = QSpacerItem(150,31,QSizePolicy.Expanding,QSizePolicy.Minimum)
        buttonsFrameLayout.addItem(spacer1)

        self.okButton = QPushButton(self.buttonsFrame,"okButton")
        self.okButton.setEnabled(0)
        buttonsFrameLayout.addWidget(self.okButton)

        self.cancelButton = QPushButton(self.buttonsFrame,"cancelButton")
        buttonsFrameLayout.addWidget(self.cancelButton)
        LayerLinkBaseLayout.addWidget(self.buttonsFrame)

        self.languageChange()

        self.resize(QSize(600,426).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.okButton,SIGNAL("clicked()"),self.accept)
        self.connect(self.cancelButton,SIGNAL("clicked()"),self.reject)
        self.connect(self.linkButton,SIGNAL("clicked()"),self.my_link)


    def languageChange(self):
        self.setCaption(self.__tr("Geometry Layer Linking"))
        QToolTip.add(self,self.__tr("Inter-Layer link generation dialog"))
        QWhatsThis.add(self,self.__tr("<h1>Inter-Layer Link Generation Dialog</h1>\n"
"Dialog to generate linkatoms between layers of the loaded geometry."))
        self.groupBox1.setTitle(self.__tr("linking parameters"))
        QToolTip.add(self.frame5,self.__tr("scaling factor to apply to link-atom equilibrium distance"))
        self.textLabel1.setText(self.__tr("link-atom distance factor"))
        self.LADFlineEdit.setText(self.__tr("1.000"))
        self.LADFlineEdit.setInputMask(self.__tr("9.000; "))
        self.textLabel2.setText(self.__tr("QMZ layer"))
        self.textLabel2_2.setText(self.__tr("external charges layer"))
        self.groupBox2.setTitle(self.__tr("linking results"))
        self.linkButton.setText(self.__tr("link"))
        self.okButton.setText(self.__tr("OK"))
        self.cancelButton.setText(self.__tr("cancel"))


    def setGeometry(self,a0):
        #make a deep copy of the original geometry object to avoid overwriting data
        	self.geo=copy.deepcopy(a0)
        # compile list of layers in Geometry
        	layers={}
        	for i in self.geo.LayerDict.keys():
        		layers[i]=self.geo.LayerDict[i].Name
        # set combobox selections
        	for i in layers.keys():
        		self.QMZComboBox.insertItem(layers[i],i)
        		self.PCHRComboBox.insertItem(layers[i],i)
        # set defaults for combocox selections: Layer QMZ for QM-zone, PCHR for external charges, if not present use default
        	if self.geo.layerbyname("QMZ")!=None:
        		self.QMZComboBox.setCurrentItem(self.geo.layerbyname("QMZ"))
        	else:
        		self.QMZComboBox.setCurrentItem(0)
        	if self.geo.layerbyname("PCHR")!=None:
        		self.PCHRComboBox.setCurrentItem(self.geo.layerbyname("PCHR"))
        	else:
        		self.PCHRComboBox.setCurrentItem(0)
        # for some reason, setting the default text of the LADF lineedit doesn't work, so set it here manually
        	self.LADFlineEdit.setText("1.00")
        		

    def my_link(self):
        	# all the outer stuff which should be common to all linking operations
        	# the real linking is done by self.linkAction called from here
        	# user input check: if both layers are the same, pop-up an error message and return without action
        	if self.QMZComboBox.currentItem()==self.PCHRComboBox.currentItem():
        		QMessageBox.warning(self,"Layer selection error","QMZ and external charges layers must be different. If your geometry has only one layer, create a new one and move some atoms into it.",QMessageBox.Ok ,QMessageBox.NoButton ,QMessageBox.NoButton)
        		return
        	#call the actual linking action
        	results=self.linkAction()
        	# put some text in results pane
        	self.resultsTextEdit.setText(utils.dictionaryPrettyPrint(results))
        	# set new button availabilities
        	self.okButton.setEnabled(True)
        	self.linkButton.setEnabled(False)
        

    def __tr(self,s,c = None):
        return qApp.translate("layerLinkBase",s,c)
