# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/knaup/src/comatsci/geostatspack/linkLayersBCTC.ui'
#
# Created: Di Sep 2 16:37:13 2008
#      by: The PyQt User Interface Compiler (pyuic) 3.17.3
#
# WARNING! All changes made in this file will be lost!


import sys
from qt import *
from qttable import QTable

from comatsci import utils
import copy


class layerLinkBCTC(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QDialog.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("LayerLinkBCTC")

        self.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Minimum,0,0,self.sizePolicy().hasHeightForWidth()))
        self.setMinimumSize(QSize(600,600))
        self.setBaseSize(QSize(600,600))
        self.setModal(1)

        LayerLinkBCTCLayout = QVBoxLayout(self,6,6,"LayerLinkBCTCLayout")

        self.parametersGroupBox = QGroupBox(self,"parametersGroupBox")
        self.parametersGroupBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Preferred,0,0,self.parametersGroupBox.sizePolicy().hasHeightForWidth()))
        self.parametersGroupBox.setMinimumSize(QSize(590,150))
        self.parametersGroupBox.setColumnLayout(0,Qt.Vertical)
        self.parametersGroupBox.layout().setSpacing(6)
        self.parametersGroupBox.layout().setMargin(6)
        parametersGroupBoxLayout = QVBoxLayout(self.parametersGroupBox.layout())
        parametersGroupBoxLayout.setAlignment(Qt.AlignTop)

        self.frame5 = QFrame(self.parametersGroupBox,"frame5")
        self.frame5.setMaximumSize(QSize(580,32767))
        self.frame5.setFrameShape(QFrame.NoFrame)
        self.frame5.setFrameShadow(QFrame.Raised)
        frame5Layout = QHBoxLayout(self.frame5,6,6,"frame5Layout")

        self.textLabel1 = QLabel(self.frame5,"textLabel1")
        frame5Layout.addWidget(self.textLabel1)

        self.LADFlineEdit = QLineEdit(self.frame5,"LADFlineEdit")
        frame5Layout.addWidget(self.LADFlineEdit)
        spacer2 = QSpacerItem(350,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
        frame5Layout.addItem(spacer2)
        parametersGroupBoxLayout.addWidget(self.frame5)

        self.frame6 = QFrame(self.parametersGroupBox,"frame6")
        self.frame6.setMaximumSize(QSize(580,32767))
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
        parametersGroupBoxLayout.addWidget(self.frame6)

        self.BCTCFrame = QFrame(self.parametersGroupBox,"BCTCFrame")
        self.BCTCFrame.setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Preferred,0,0,self.BCTCFrame.sizePolicy().hasHeightForWidth()))
        self.BCTCFrame.setMinimumSize(QSize(560,0))
        self.BCTCFrame.setMaximumSize(QSize(580,32767))
        self.BCTCFrame.setFrameShape(QFrame.StyledPanel)
        self.BCTCFrame.setFrameShadow(QFrame.Raised)
        BCTCFrameLayout = QGridLayout(self.BCTCFrame,1,1,6,6,"BCTCFrameLayout")

        self.BCTCCalculateButton = QPushButton(self.BCTCFrame,"BCTCCalculateButton")

        BCTCFrameLayout.addWidget(self.BCTCCalculateButton,2,1)

        self.BCTCFileButton = QPushButton(self.BCTCFrame,"BCTCFileButton")

        BCTCFrameLayout.addWidget(self.BCTCFileButton,2,0)

        self.BCTCTable = QTable(self.BCTCFrame,"BCTCTable")
        self.BCTCTable.setNumCols(self.BCTCTable.numCols() + 1)
        self.BCTCTable.horizontalHeader().setLabel(self.BCTCTable.numCols() - 1,self.__tr("BCT"))
        self.BCTCTable.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.MinimumExpanding,0,0,self.BCTCTable.sizePolicy().hasHeightForWidth()))
        self.BCTCTable.setMinimumSize(QSize(250,100))
        self.BCTCTable.setMaximumSize(QSize(260,32767))
        self.BCTCTable.setNumRows(0)
        self.BCTCTable.setNumCols(1)
        self.BCTCTable.setSorting(1)

        BCTCFrameLayout.addMultiCellWidget(self.BCTCTable,0,2,2,2)

        self.textLabel1_2 = QLabel(self.BCTCFrame,"textLabel1_2")

        BCTCFrameLayout.addWidget(self.textLabel1_2,1,0)

        self.neutralizeClusterCheckBox = QCheckBox(self.BCTCFrame,"neutralizeClusterCheckBox")
        self.neutralizeClusterCheckBox.setChecked(1)

        BCTCFrameLayout.addMultiCellWidget(self.neutralizeClusterCheckBox,0,0,0,1)
        parametersGroupBoxLayout.addWidget(self.BCTCFrame)
        LayerLinkBCTCLayout.addWidget(self.parametersGroupBox)

        self.groupBox2 = QGroupBox(self,"groupBox2")
        self.groupBox2.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Expanding,0,0,self.groupBox2.sizePolicy().hasHeightForWidth()))
        self.groupBox2.setMinimumSize(QSize(590,200))
        self.groupBox2.setColumnLayout(0,Qt.Vertical)
        self.groupBox2.layout().setSpacing(6)
        self.groupBox2.layout().setMargin(6)
        groupBox2Layout = QVBoxLayout(self.groupBox2.layout())
        groupBox2Layout.setAlignment(Qt.AlignTop)

        self.resultsTextEdit = QTextEdit(self.groupBox2,"resultsTextEdit")
        self.resultsTextEdit.setEnabled(1)
        self.resultsTextEdit.setTextFormat(QTextEdit.PlainText)
        self.resultsTextEdit.setWordWrap(QTextEdit.NoWrap)
        self.resultsTextEdit.setReadOnly(1)
        self.resultsTextEdit.setUndoRedoEnabled(0)
        groupBox2Layout.addWidget(self.resultsTextEdit)
        LayerLinkBCTCLayout.addWidget(self.groupBox2)

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
        LayerLinkBCTCLayout.addWidget(self.buttonsFrame)

        self.languageChange()

        self.resize(QSize(600,677).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.okButton,SIGNAL("clicked()"),self.accept)
        self.connect(self.cancelButton,SIGNAL("clicked()"),self.reject)
        self.connect(self.linkButton,SIGNAL("clicked()"),self.my_link)
        self.connect(self.BCTCCalculateButton,SIGNAL("clicked()"),self.BCTCCoefficientsFromCalc)
        self.connect(self.BCTCFileButton,SIGNAL("clicked()"),self.BCTCCoefficientsFromFile)


    def languageChange(self):
        self.setCaption(self.__tr("Geometry Layer Linking"))
        QToolTip.add(self,self.__tr("Inter-Layer link generation dialog"))
        QWhatsThis.add(self,self.__tr("<h1>Inter-Layer Link Generation Dialog</h1>\n"
"Dialog to generate linkatoms between layers of the loaded geometry."))
        self.parametersGroupBox.setTitle(self.__tr("linking parameters"))
        QToolTip.add(self.frame5,self.__tr("scaling factor to apply to link-atom equilibrium distance"))
        self.textLabel1.setText(self.__tr("link-atom distance factor"))
        self.LADFlineEdit.setText(self.__tr("1.00"))
        self.LADFlineEdit.setInputMask(self.__tr("9.000; "))
        self.textLabel2.setText(self.__tr("QMZ layer"))
        self.textLabel2_2.setText(self.__tr("external charges layer"))
        self.BCTCCalculateButton.setText(self.__tr("calculate"))
        self.BCTCFileButton.setText(self.__tr("read file..."))
        self.BCTCTable.horizontalHeader().setLabel(0,self.__tr("BCT"))
        QToolTip.add(self.BCTCTable,self.__tr("Bonc Charge Transfer Coefficients"))
        self.textLabel1_2.setText(self.__tr("BCT coefficients"))
        self.neutralizeClusterCheckBox.setText(self.__tr("neutralize cluster"))
        QToolTip.add(self.neutralizeClusterCheckBox,self.__tr("neutralize the resildual cluster charge after BCT compensation by HCS"))
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
        #generate initial BCT coefficients dictionary with zero coefficients and populate BCTCTable
        	(atomsymbols,dummy)=self.geo.getatomsymlistdict()
        	atomsymbols.sort()
        	BCTCoefficients={}
        	for i in range(len(atomsymbols)):
        		for j in range(i,len(atomsymbols)):
        			BCTCoefficients[(atomsymbols[i],atomsymbols[j])]=0.0
        	self.setBCTCTable(BCTCoefficients)
        		

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
        

    def setBCTCTable(self,a0):
        	#clear the BCTCTable and populate it with the given BCTC coefficients
        	# clear the table
        	self.BCTCTable.removeRows(range(self.BCTCTable.numRows()))
        	# add the new rows
        	combinations=a0.keys()
        	combinations.sort()
        	self.BCTCTable.insertRows(0,len(combinations))
        	rowlabels=QStringList()
        	for row in range(len(combinations)):
        		self.BCTCTable.setText(row,0,str(a0[combinations[row]]))
        		rowlabels.append("%d-%d"%combinations[row])
        	# apply the row labels and make column fill the available space
        	self.BCTCTable.setRowLabels(rowlabels)
        	self.BCTCTable.adjustColumn(0)
        	# store the BCTC coeffient key<->BCTCTable row mapping
        	self.BCTCTableKeys=combinations
        

    def BCTCCoefficientsFromTable(self):
        	BCTCoefficients={}
        	# get the BCTC values from the table and return them in comatsci-friendly format
        	for i in range(len(self.BCTCTableKeys)):
        		try:
        			BCTCoefficients[self.BCTCTableKeys[i]]=float(str(self.BCTCTable.text(i,0)))
        		except:
        			QMessageBox.warning(self,"BCTC Coefficient parse error","The BCTC  coefficient specified in line %s of the coefficients table could not be parsed. Aborting."%(self.BCTCTableKeys[i],),QMessageBox.Cancel ,QMessageBox.NoButton ,QMessageBox.NoButton)
        			raise
        	# return coefficients dictionary
        	return BCTCoefficients
        

    def BCTCCoefficientsFromFile(self):
        	#read the BCTCCoefficients from a file and store them in BCTCTable
        	# pop-up modal file selecttion dialog
        	filename=str(QFileDialog.getOpenFileName(".","Known formats (*.bctc *.bctc.gz)", self))
        	if (len(filename)==0):
        		return
        	# pull input file to memory
        	try:
        		infile=utils.compressedopen(filename)
        		inlines=list(infile)
        		infile.close()
        	except:
        		QMessageBox.warning(self,"BCTC Coefficient read error","The specified BCTC  coefficient file '%s' could not be read. Aborting."%(filename,),QMessageBox.Cancel ,QMessageBox.NoButton ,QMessageBox.NoButton)
        		raise
        	# initialize BCTC coefficients dictionary
        	BCTCCoefficients={}
        	# parse input file lines. ignore marked comment lines and give warnings to the user if line parsing fails
        	for i in range(len(inlines)):
        		if not inlines[i].strip()[0] in [';','#']:
        			try:
        				parts=inlines[i].split()
        				Z1=int(parts[0])
        				Z2=int(parts[1])
        				C=float(parts[2])
        			except:
        				QMessageBox.warning(self,"BCTC Coefficient read error","Fatal error parsing line %d of file '%s'. Expecting lines of two integers and one float. Aborting."%(i+1,filename),QMessageBox.Cancel ,QMessageBox.NoButton ,QMessageBox.NoButton)
        				raise
        			BCTCCoefficients[(Z1,Z2)]=C
        	# store the Coefficients into the dialog's table
        	self.setBCTCTable(BCTCCoefficients)
        

    def BCTCCoefficientsFromCalc(self):
        	# calculate BCTC coefficients from geometry and apply to the coefficients table
        	BCTCoefficients=self.geo.getElementElementChargeTransfers()[0]
        	self.setBCTCTable(BCTCoefficients)
        

    def linkAction(self):
        		# get the link-atom distance factor
        		ladf=float(self.LADFlineEdit.text())
        		# get the BCTC coefficients
        		coefficients=self.BCTCCoefficientsFromTable()
        		# create a progress dialog and execute the linking operation
        		progressWindow=QProgressDialog()
        		progressWindow.setLabelText("generating BCTC atoms")
        		progressWindow.show()
        		(self.embeddedGeometry,results)=self.geo.layersubgeometry(self.QMZComboBox.currentItem()).BCTCLinkedGeometry(self.geo.layersubgeometry(self.PCHRComboBox.currentItem()),progressWindow.setTotalSteps,progressWindow.setProgress,distscale=ladf,chargeTransfers=coefficients,neutralize=self.neutralizeClusterCheckBox.isChecked())
        		progressWindow.hide()
        		return results
        

    def __tr(self,s,c = None):
        return qApp.translate("layerLinkBCTC",s,c)

if __name__ == "__main__":
    a = QApplication(sys.argv)
    QObject.connect(a,SIGNAL("lastWindowClosed()"),a,SLOT("quit()"))
    w = layerLinkBCTC()
    a.setMainWidget(w)
    w.show()
    a.exec_loop()
