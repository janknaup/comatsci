# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'rdfWizard.ui'
#
# Created: Do Jul 16 16:28:05 2009
#      by: The PyQt User Interface Compiler (pyuic) 3.17.6
#
# WARNING! All changes made in this file will be lost!


from qt import *
##############################################################################
# rdfWizard.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup , Knaup@bccms.uni-bremen.de
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


class rdfWizard(QWizard):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QWizard.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("rdfWizard")



        self.page = QWidget(self,"page")
        pageLayout = QGridLayout(self.page,1,1,6,6,"pageLayout")

        self.textLabel1 = QLabel(self.page,"textLabel1")

        pageLayout.addWidget(self.textLabel1,1,0)

        self.textLabel1_3 = QLabel(self.page,"textLabel1_3")

        pageLayout.addWidget(self.textLabel1_3,2,2)

        self.stepWidthLCDNumber = QLCDNumber(self.page,"stepWidthLCDNumber")
        self.stepWidthLCDNumber.setFrameShadow(QLCDNumber.Sunken)
        self.stepWidthLCDNumber.setNumDigits(4)
        self.stepWidthLCDNumber.setSegmentStyle(QLCDNumber.Filled)

        pageLayout.addWidget(self.stepWidthLCDNumber,2,1)

        self.stepWidthSlider = QSlider(self.page,"stepWidthSlider")
        self.stepWidthSlider.setMaxValue(200)
        self.stepWidthSlider.setValue(20)
        self.stepWidthSlider.setOrientation(QSlider.Horizontal)
        self.stepWidthSlider.setTickmarks(QSlider.Above)
        self.stepWidthSlider.setTickInterval(25)

        pageLayout.addWidget(self.stepWidthSlider,2,0)

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

        self.textLabel1_2 = QLabel(self.WizardPage_2,"textLabel1_2")

        WizardPageLayout_2.addWidget(self.textLabel1_2,1,0)

        self.textLabel2_2 = QLabel(self.WizardPage_2,"textLabel2_2")

        WizardPageLayout_2.addWidget(self.textLabel2_2,3,0)

        self.textLabel7 = QLabel(self.WizardPage_2,"textLabel7")

        WizardPageLayout_2.addWidget(self.textLabel7,0,0)

        self.postProcessProgressBar = QProgressBar(self.WizardPage_2,"postProcessProgressBar")

        WizardPageLayout_2.addWidget(self.postProcessProgressBar,4,0)

        self.binningProgressBar = QProgressBar(self.WizardPage_2,"binningProgressBar")

        WizardPageLayout_2.addWidget(self.binningProgressBar,2,0)
        self.addPage(self.WizardPage_2,QString(""))

        self.languageChange()

        self.resize(QSize(390,358).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.callGnuplotCheckBox,SIGNAL("toggled(bool)"),self.setSavePlotState)
        self.connect(self.stepWidthSlider,SIGNAL("valueChanged(int)"),self.stepWidthChanged)

        self.setTabOrder(self.dataFileEdit,self.callGnuplotCheckBox)
        self.setTabOrder(self.callGnuplotCheckBox,self.savePlotCheckBox)
        self.setTabOrder(self.savePlotCheckBox,self.plotFileEdit)

        self.init()


    def languageChange(self):
        self.setCaption(self.__tr("radial distribution function","Bin atomic charges into historams per element"))
        QToolTip.add(self,self.__tr("Wizard to calculate the radial distribution function"))
        QWhatsThis.add(self,self.__tr("<h1>Radial distribution function (rdf) wizard</h1>\n"
"<p>The rdf wizard helps you to calculate the radial distribution function of the given Geometry. For supercells, the cell volume is used for normalization, otherwise a sphere thightly circumcising the molecule. <b>Normalization breaks for slab and wire models, due to the vacuum in the supercell.<b></p>\n"
"<p>The rdf is saved into a fixed-column ascii datafile.</p>\n"
"<p>If the Python Gnuplot package and gnuplot are available, the rdf can be directly plotted and the plot can be saved.</p>"))
        self.textLabel1.setText(self.__tr("<b>Spherical shell thickness</b>"))
        self.textLabel1_3.setText(self.__tr("<font size=\"+3\">Bohr</font>"))
        self.textLabel2.setText(self.__tr("<h1>Stepwidth</h1>\n"
"<p>The radial distribution function is calculated by binning the interatomic distances into spherical shells of equal thickness and dividing the counts by sphere shell volume.</p>\n"
"<p>Select the rdf stepwidth in Bohr radii.</p>"))
        self.setTitle(self.page,self.__tr("Binning"))
        self.textLabel4.setText(self.__tr("Data file name"))
        self.dataFileEdit.setText(self.__tr("rdf.dat"))
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
        self.textLabel1_2.setText(self.__tr("<b>Binning</b>"))
        self.textLabel2_2.setText(self.__tr("<b>Postprocessing</b>"))
        self.textLabel7.setText(self.__tr("<h1>Finishing</h1>\n"
"<p>Binning scales like N<sup>2</sup> with atom count, the Posprocessing scales linear.</p>"))
        self.setTitle(self.WizardPage_2,self.__tr("Finish"))


    def setSavePlotState(self,a0):
        	self.savePlotCheckBox.setEnabled(a0)
        	self.plotFileEdit.setEnabled(a0)
        	self.plotFilenameTextLabel.setEnabled(a0)
        	self.fileFormatButtonGroup.setEnabled(a0)
        

    def accept(self):
        	self.rdfdata=self.geo.rdf(self.stepwidth, 
        		self.postProcessProgressBar.setTotalSteps, self.postProcessProgressBar.setProgress,
        		self.binningProgressBar.setTotalSteps, self.binningProgressBar.setProgress)
        	datafile=open(str(self.dataFileEdit.text()),'w')
        	print >> datafile, "# r[a.u.]  rdf"
        	for j in range(len(self.rdfdata[0])):
        		print >> datafile, "%e   %e" % (self.rdfdata[0][j],self.rdfdata[1][j])
        	datafile.close()
        	# calculate elemental RDFs and write datafiles
        	elementrdfs=self.geo.elementRDFs(self.stepwidth, self.postProcessProgressBar.setTotalSteps, self.postProcessProgressBar.setProgress,
        		self.binningProgressBar.setTotalSteps, self.binningProgressBar.setProgress)
        	for i in elementrdfs.keys():
        		elementDataFile=open(str(self.dataFileEdit.text())[0:-4]+"-"+self.geo.PTE[i]+str(self.dataFileEdit.text())[-4:],"w")
        		print >> elementDataFile, "# r[a.u.] rdf_Z"
        		for j in range(len(elementrdfs[i][0])):
        			print >> elementDataFile, "%e   %e" % (elementrdfs[i][0][j],elementrdfs[i][1][j])
        		elementDataFile.close()
        	if self.callGnuplotCheckBox.isChecked():
        		self.plotHistogram()
        	QDialog.accept(self)
        

    def plotHistogram(self):
        	rdfplot=self.Gnuplot.Data(self.rdfdata[0],self.rdfdata[1],title="total rdf",inline=1)
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
        		self.plot("set title 'total radial distribution function'")
        		self.plot("set xlabel 'r[Bohr]'")
        		self.plot("set ylabel 'rdf'")
        	#initialize stepwidth input slider and LCD display
        	self.stepWidthSlider.setValue(20)
        	self.stepWidthLCDNumber.display(0.20)
        	self.stepwidth=0.20
        

    def stepWidthChanged(self,a0):
        	self.stepwidth=float(a0)/100.0
        	self.stepWidthLCDNumber.display(self.stepwidth)
        

    def __tr(self,s,c = None):
        return qApp.translate("rdfWizard",s,c)
