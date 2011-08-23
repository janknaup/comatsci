# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'anglesWizard.ui'
#
# Created: Wed Aug 17 14:31:16 2011
#      by: The PyQt User Interface Compiler (pyuic) 3.18.1
#
# WARNING! All changes made in this file will be lost!


from qt import *
##############################################################################
# anglesWizard.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup , Knaup@bccms.uni-bremen.de
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


class anglesWizard(QWizard):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QWizard.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("angleWizard")



        self.page = QWidget(self,"page")
        pageLayout = QGridLayout(self.page,1,1,6,6,"pageLayout")

        self.textLabel1 = QLabel(self.page,"textLabel1")

        pageLayout.addWidget(self.textLabel1,1,0)

        self.textLabel1_3 = QLabel(self.page,"textLabel1_3")

        pageLayout.addWidget(self.textLabel1_3,2,2)

        self.binsLCDNumber = QLCDNumber(self.page,"binsLCDNumber")
        self.binsLCDNumber.setFrameShadow(QLCDNumber.Sunken)
        self.binsLCDNumber.setNumDigits(4)
        self.binsLCDNumber.setSegmentStyle(QLCDNumber.Filled)

        pageLayout.addWidget(self.binsLCDNumber,2,1)

        self.binsSlider = QSlider(self.page,"binsSlider")
        self.binsSlider.setMaxValue(200)
        self.binsSlider.setValue(20)
        self.binsSlider.setOrientation(QSlider.Horizontal)
        self.binsSlider.setTickmarks(QSlider.Above)
        self.binsSlider.setTickInterval(25)

        pageLayout.addWidget(self.binsSlider,2,0)

        self.textLabel2 = QLabel(self.page,"textLabel2")

        pageLayout.addMultiCellWidget(self.textLabel2,0,0,0,2)
        self.addPage(self.page,QString(""))

        self.WizardPage = QWidget(self,"WizardPage")
        WizardPageLayout = QGridLayout(self.WizardPage,1,1,6,6,"WizardPageLayout")

        self.textLabel4 = QLabel(self.WizardPage,"textLabel4")

        WizardPageLayout.addWidget(self.textLabel4,1,0)

        self.dataFileEdit = QLineEdit(self.WizardPage,"dataFileEdit")

        WizardPageLayout.addWidget(self.dataFileEdit,1,1)

        self.textLabel3 = QLabel(self.WizardPage,"textLabel3")

        WizardPageLayout.addMultiCellWidget(self.textLabel3,0,0,0,1)

        self.frame3 = QFrame(self.WizardPage,"frame3")
        self.frame3.setFrameShape(QFrame.StyledPanel)
        self.frame3.setFrameShadow(QFrame.Raised)
        frame3Layout = QGridLayout(self.frame3,1,1,6,6,"frame3Layout")

        self.callGnuplotCheckBox = QCheckBox(self.frame3,"callGnuplotCheckBox")
        self.callGnuplotCheckBox.setTristate(0)

        frame3Layout.addMultiCellWidget(self.callGnuplotCheckBox,0,1,0,1)

        self.plotFilenameTextLabel = QLabel(self.frame3,"plotFilenameTextLabel")
        self.plotFilenameTextLabel.setEnabled(0)

        frame3Layout.addWidget(self.plotFilenameTextLabel,2,1)

        self.plotFileEdit = QLineEdit(self.frame3,"plotFileEdit")
        self.plotFileEdit.setEnabled(0)

        frame3Layout.addWidget(self.plotFileEdit,3,1)

        self.savePlotCheckBox = QCheckBox(self.frame3,"savePlotCheckBox")
        self.savePlotCheckBox.setEnabled(0)

        frame3Layout.addWidget(self.savePlotCheckBox,1,1)

        self.fileFormatButtonGroup = QButtonGroup(self.frame3,"fileFormatButtonGroup")
        self.fileFormatButtonGroup.setEnabled(0)
        self.fileFormatButtonGroup.setFlat(0)

        self.emfRadioButton = QRadioButton(self.fileFormatButtonGroup,"emfRadioButton")
        self.emfRadioButton.setGeometry(QRect(10,40,60,21))

        self.pngRadioButton = QRadioButton(self.fileFormatButtonGroup,"pngRadioButton")
        self.pngRadioButton.setGeometry(QRect(80,20,60,21))

        self.figRadioButton = QRadioButton(self.fileFormatButtonGroup,"figRadioButton")
        self.figRadioButton.setGeometry(QRect(80,40,60,21))

        self.epsRadioButton = QRadioButton(self.fileFormatButtonGroup,"epsRadioButton")
        self.epsRadioButton.setGeometry(QRect(10,20,60,21))
        self.epsRadioButton.setChecked(1)

        frame3Layout.addMultiCellWidget(self.fileFormatButtonGroup,1,3,2,2)

        WizardPageLayout.addMultiCellWidget(self.frame3,2,2,0,1)
        self.addPage(self.WizardPage,QString(""))

        self.WizardPage_2 = QWidget(self,"WizardPage_2")
        WizardPageLayout_2 = QGridLayout(self.WizardPage_2,1,1,6,6,"WizardPageLayout_2")

        self.textLabel7 = QLabel(self.WizardPage_2,"textLabel7")

        WizardPageLayout_2.addWidget(self.textLabel7,0,0)
        self.addPage(self.WizardPage_2,QString(""))

        self.languageChange()

        self.resize(QSize(390,358).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.callGnuplotCheckBox,SIGNAL("toggled(bool)"),self.setSavePlotState)
        self.connect(self.binsSlider,SIGNAL("valueChanged(int)"),self.binsChanged)

        self.setTabOrder(self.dataFileEdit,self.callGnuplotCheckBox)
        self.setTabOrder(self.callGnuplotCheckBox,self.savePlotCheckBox)
        self.setTabOrder(self.savePlotCheckBox,self.plotFileEdit)

        self.init()


    def languageChange(self):
        self.setCaption(self.__tr("bond angle histograms","Bin atomic charges into historams per element"))
        QToolTip.add(self,self.__tr("Wizard to calculate the radial distribution function"))
        QWhatsThis.add(self,self.__tr("<h1>Radial distribution function (rdf) wizard</h1>\n"
"<p>The rdf wizard helps you to calculate the radial distribution function of the given Geometry. For supercells, the cell volume is used for normalization, otherwise a sphere thightly circumcising the molecule. <b>Normalization breaks for slab and wire models, due to the vacuum in the supercell.<b></p>\n"
"<p>The rdf is saved into a fixed-column ascii datafile.</p>\n"
"<p>If the Python Gnuplot package and gnuplot are available, the rdf can be directly plotted and the plot can be saved.</p>"))
        self.textLabel1.setText(self.__tr("<b>Bin Count</b>"))
        self.textLabel1_3.setText(self.__tr("<font size=\"+3\">Bins</font>"))
        self.textLabel2.setText(self.__tr("<h1>Stepwidth</h1>\n"
"<p>Histograms of bond angle distributions by central atom elements.</p>\n"
"<p>Select the number of Bins.</p>"))
        self.setTitle(self.page,self.__tr("Binning"))
        self.textLabel4.setText(self.__tr("Data file name"))
        self.dataFileEdit.setText(self.__tr("anglehist.dat"))
        self.textLabel3.setText(self.__tr("<h1>Output options</h1>\n"
"<p>Select the data filename, decide wether to call gnuplot and to directly save the plot. Choose the plot filename and -type if applicable.</p>"))
        self.callGnuplotCheckBox.setText(self.__tr("call gnuplot"))
        self.plotFilenameTextLabel.setText(self.__tr("filename prefix"))
        self.plotFileEdit.setText(self.__tr("rdf"))
        self.savePlotCheckBox.setText(self.__tr("save plot"))
        self.fileFormatButtonGroup.setTitle(self.__tr("graphics file format"))
        self.emfRadioButton.setText(self.__tr("emf"))
        self.pngRadioButton.setText(self.__tr("png"))
        self.figRadioButton.setText(self.__tr("fig"))
        self.epsRadioButton.setText(self.__tr("eps"))
        self.setTitle(self.WizardPage,self.__tr("Output"))
        self.textLabel7.setText(self.__tr("<h1>Finishing</h1>"))
        self.setTitle(self.WizardPage_2,self.__tr("Finish"))


    def setSavePlotState(self,a0):
        	self.savePlotCheckBox.setEnabled(a0)
        	self.plotFileEdit.setEnabled(a0)
        	self.plotFilenameTextLabel.setEnabled(a0)
        	self.fileFormatButtonGroup.setEnabled(a0)
        

    def accept(self):
        	self.totalHistogram=self.geo.getBondAngleHistogram(self.bins)
        	datafile=open(str(self.dataFileEdit.text()),'w')
        	print >> datafile, "# angle  total count"
        	for j in range(len(self.totalHistogram[0])):
        		print >> datafile, "%e   %e" % (self.totalHistogram[1][j],self.totalHistogram[0][j])
        	datafile.close()
        	# calculate elemental RDFs and write datafiles
        	elementHistograms=self.geo.getBondAnglesByElementsHistograms(self.bins)
        	for i in elementHistograms.keys():
        		elementDataFile=open(str(self.dataFileEdit.text())[0:-4]+"-"+self.geo.PTE[i]+str(self.dataFileEdit.text())[-4:],"w")
        		print >> elementDataFile, "# angle count"
        		for j in range(len(elementHistograms[i][0])):
        			print >> elementDataFile, "%e   %e" % (elementHistograms[i][1][j],elementHistograms[i][0][j])
        		elementDataFile.close()
        	if self.callGnuplotCheckBox.isChecked():
        		self.plotHistogram()
        	QDialog.accept(self)
        

    def plotHistogram(self):
        	rdfplot=self.Gnuplot.Data(self.totalHistogram[0],self.totHistograms[1],title="total angular distribution",inline=1)
        	self.plot.plot(rdfplot)
        	if self.savePlotCheckBox.isChecked():
        		savefilename=str(self.plotFileEdit.text())
        		if self.emfRadioButton.isChecked():
        			self.plot("set term emf color solid")
        			if savefilename[-4:]!=".emf":
        				savefilename+=".emf"
        		elif self.epsRadioButton.isChecked():
        			self.plot("set term post eps col sol enh")
        			if savefilename[-4:]!=".eps":
        				savefilename+=".eps"
        		elif self.pngRadioButton.isChecked():
        			self.plot("set term png enhanced")
        			if savefilename[-4:]!=".png":
        				savefilename+=".png"
        		elif self.figRadioButton.isChecked():
        			self.plot("set term fig col sol")
        			if savefilename[-4:]!=".fig":
        				savefilename+=".fig"
        		else:
        			pass
        		self.plot("set out '%s'" % (savefilename))
        		self.plot.plot(rdfplot)
        		self.plot("set out")
        

    def setGeometry(self,a0):
        	self.geo=a0
        	self.setFinishEnabled(self.WizardPage_2,True)
        

    def init(self):
        	"""Initialize Charge histogram wizard
        	check for Gnuplot package and init that if available."""
        	try:
        		self.Gnuplot= __import__("Gnuplot")
        	except:
        		self.callGnuplotCheckBox.setEnabled(False)
        	else:
        		self.callGnuplotCheckBox.setEnabled(True)
        		self.plot=self.Gnuplot.Gnuplot(debug=0,persist=1)
        		self.plot("set data style lines")
        		self.plot("set title 'total angular distribution'")
        		self.plot("set xlabel 'angle'")
        		self.plot("set ylabel 'count'")
        	#initialize bins input slider and LCD display
        	self.binsSlider.setValue(36)
        	self.binsLCDNumber.display(36)
        	self.bins=36
        

    def binsChanged(self,a0):
        	self.bins=int(a0)
        	self.binsLCDNumber.display(self.bins)
        

    def __tr(self,s,c = None):
        return qApp.translate("anglesWizard",s,c)
