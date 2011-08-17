# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'elementChargeHistograms.ui'
#
# Created: Wed Aug 17 14:31:16 2011
#      by: The PyQt User Interface Compiler (pyuic) 3.18.1
#
# WARNING! All changes made in this file will be lost!


from qt import *
##############################################################################
# elementChargeHistograms.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup , Knaup@bccms.uni-bremen.de
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


class elementChargesHistograms(QWizard):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QWizard.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("Form")



        self.page = QWidget(self,"page")

        self.textLabel2 = QLabel(self.page,"textLabel2")
        self.textLabel2.setGeometry(QRect(10,0,370,130))

        self.textLabel1 = QLabel(self.page,"textLabel1")
        self.textLabel1.setGeometry(QRect(170,140,210,20))

        self.binsSpinBox = QSpinBox(self.page,"binsSpinBox")
        self.binsSpinBox.setGeometry(QRect(0,140,151,20))
        self.binsSpinBox.setMaxValue(1000)
        self.binsSpinBox.setMinValue(2)
        self.binsSpinBox.setValue(10)
        self.addPage(self.page,QString(""))

        self.WizardPage = QWidget(self,"WizardPage")

        self.textLabel4 = QLabel(self.WizardPage,"textLabel4")
        self.textLabel4.setGeometry(QRect(10,110,100,20))

        self.dataFileEdit = QLineEdit(self.WizardPage,"dataFileEdit")
        self.dataFileEdit.setGeometry(QRect(120,110,240,21))

        self.textLabel3 = QLabel(self.WizardPage,"textLabel3")
        self.textLabel3.setGeometry(QRect(10,0,351,100))

        self.frame3_2 = QFrame(self.WizardPage,"frame3_2")
        self.frame3_2.setGeometry(QRect(10,140,360,100))
        self.frame3_2.setFrameShape(QFrame.StyledPanel)
        self.frame3_2.setFrameShadow(QFrame.Raised)

        self.callGnuplotCheckBox = QCheckBox(self.frame3_2,"callGnuplotCheckBox")
        self.callGnuplotCheckBox.setGeometry(QRect(6,3,100,20))
        self.callGnuplotCheckBox.setTristate(0)

        self.savePlotCheckBox = QCheckBox(self.frame3_2,"savePlotCheckBox")
        self.savePlotCheckBox.setEnabled(0)
        self.savePlotCheckBox.setGeometry(QRect(20,20,151,21))

        self.plotFilenameTextLabel = QLabel(self.frame3_2,"plotFilenameTextLabel")
        self.plotFilenameTextLabel.setEnabled(0)
        self.plotFilenameTextLabel.setGeometry(QRect(20,50,93,21))

        self.fileFormatButtonGroup = QButtonGroup(self.frame3_2,"fileFormatButtonGroup")
        self.fileFormatButtonGroup.setEnabled(0)
        self.fileFormatButtonGroup.setGeometry(QRect(190,20,150,70))
        self.fileFormatButtonGroup.setFlat(0)

        self.pngRadioButton = QRadioButton(self.fileFormatButtonGroup,"pngRadioButton")
        self.pngRadioButton.setGeometry(QRect(80,20,60,21))

        self.figRadioButton = QRadioButton(self.fileFormatButtonGroup,"figRadioButton")
        self.figRadioButton.setGeometry(QRect(80,40,60,21))

        self.epsRadioButton = QRadioButton(self.fileFormatButtonGroup,"epsRadioButton")
        self.epsRadioButton.setGeometry(QRect(10,20,60,21))
        self.epsRadioButton.setChecked(1)

        self.emfRadioButton = QRadioButton(self.fileFormatButtonGroup,"emfRadioButton")
        self.emfRadioButton.setGeometry(QRect(10,40,60,21))

        self.plotFileEdit = QLineEdit(self.frame3_2,"plotFileEdit")
        self.plotFileEdit.setEnabled(0)
        self.plotFileEdit.setGeometry(QRect(20,70,140,21))
        self.addPage(self.WizardPage,QString(""))

        self.WizardPage_2 = QWidget(self,"WizardPage_2")

        self.textLabel7 = QLabel(self.WizardPage_2,"textLabel7")
        self.textLabel7.setGeometry(QRect(10,10,270,60))

        self.textLabel1_2 = QLabel(self.WizardPage_2,"textLabel1_2")
        self.textLabel1_2.setGeometry(QRect(10,90,281,16))

        self.textLabel2_2 = QLabel(self.WizardPage_2,"textLabel2_2")
        self.textLabel2_2.setGeometry(QRect(4,148,291,20))

        self.elementsProgressBar = QProgressBar(self.WizardPage_2,"elementsProgressBar")
        self.elementsProgressBar.setGeometry(QRect(0,110,370,31))

        self.binningProgressBar = QProgressBar(self.WizardPage_2,"binningProgressBar")
        self.binningProgressBar.setGeometry(QRect(0,170,370,31))
        self.addPage(self.WizardPage_2,QString(""))

        self.languageChange()

        self.resize(QSize(390,324).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.callGnuplotCheckBox,SIGNAL("toggled(bool)"),self.setSavePlotState)

        self.setTabOrder(self.binsSpinBox,self.dataFileEdit)

        self.init()


    def languageChange(self):
        self.setCaption(self.__tr("element charge histograms","Bin atomic charges into historams per element"))
        QToolTip.add(self,self.__tr("Wizard to generate element-wise atomic charge histograms"))
        QWhatsThis.add(self,self.__tr("<h1>Element-wise atomic charge histogram wizard</h1>\n"
"<p>Generates atomic (Mulliken) charge histograms for each element in the geometry.</p>\n"
"<p>The histograms are saved into a fixed-column multi-block gnuplot compatible ascii datafile.</p>\n"
"<p>If the Python Gnuplot package and gnuplot are available, the chrge histograms can be plotted directlz and the plot can be saved.</p>"))
        self.textLabel2.setText(self.__tr("<h1>Binning</h1>\n"
"<p>The Charges are sorted into evenly spaced bins distributed between the minimum and maximum atomic charge per element.</p>\n"
"<p>Select the number of bins per element charge histogram.</p>"))
        self.textLabel1.setText(self.__tr("Bins per histogram"))
        QToolTip.add(self.binsSpinBox,self.__tr("number of bins per histogram"))
        QWhatsThis.add(self.binsSpinBox,self.__tr("The charges are binned into a range of intervals between the smallest and largest atomic charge per element. The number of intervals (=bins) per element is selected here."))
        self.setTitle(self.page,self.__tr("Binning"))
        self.textLabel4.setText(self.__tr("Data file name"))
        self.dataFileEdit.setText(self.__tr("charges.hst"))
        QToolTip.add(self.dataFileEdit,self.__tr("data filename"))
        QWhatsThis.add(self.dataFileEdit,self.__tr("Filename for the gnuplot-firendly file containing the histogram data."))
        self.textLabel3.setText(self.__tr("<h1>Output options</h1>\n"
"<p>Select the data filename, decide wether to call gnuplot and to directly save the plot. Choose the plot filename and -type if applicable.</p>"))
        QToolTip.add(self.frame3_2,self.__tr("gnuplot options"))
        QWhatsThis.add(self.frame3_2,self.__tr("If gnuplot and the gnuplot python module are present on yout system, gnuplot can be called directly to display the diagrams(s) and optionally save the plot into a graphics file."))
        self.callGnuplotCheckBox.setText(self.__tr("call gnuplot"))
        QToolTip.add(self.callGnuplotCheckBox,self.__tr("Should gnuplot be called?"))
        self.savePlotCheckBox.setText(self.__tr("save plot"))
        QToolTip.add(self.savePlotCheckBox,self.__tr("Should the plot be saved as a graphics file?"))
        self.plotFilenameTextLabel.setText(self.__tr("filename prefix"))
        self.fileFormatButtonGroup.setTitle(self.__tr("graphics file format"))
        QToolTip.add(self.fileFormatButtonGroup,self.__tr("Select the graphics file format to save the plot."))
        self.pngRadioButton.setText(self.__tr("png"))
        self.figRadioButton.setText(self.__tr("fig"))
        self.epsRadioButton.setText(self.__tr("eps"))
        self.emfRadioButton.setText(self.__tr("emf"))
        self.plotFileEdit.setText(self.__tr("charges"))
        QToolTip.add(self.plotFileEdit,self.__tr("graphics file name prefix"))
        self.setTitle(self.WizardPage,self.__tr("Output"))
        self.textLabel7.setText(self.__tr("<h1>Finishing</h1>\n"
"<p>Binning charges and writing output files.</p>"))
        self.textLabel1_2.setText(self.__tr("<b>Elements</b>"))
        self.textLabel2_2.setText(self.__tr("<b>Atom binning</b>"))
        QToolTip.add(self.elementsProgressBar,self.__tr("progress finishing the histograms per element"))
        QToolTip.add(self.binningProgressBar,self.__tr("progress of binning one element's atomic charges"))
        self.setTitle(self.WizardPage_2,self.__tr("Finish"))


    def setSavePlotState(self,a0):
        	self.savePlotCheckBox.setEnabled(a0)
        	self.plotFileEdit.setEnabled(a0)
        	self.plotFilenameTextLabel.setEnabled(a0)
        	self.fileFormatButtonGroup.setEnabled(a0)
        

    def accept(self):
        	self.chrhist=self.geo.elem_charges_hist(int(str(self.binsSpinBox.value())), 
        		self.elementsProgressBar.setTotalSteps, self.elementsProgressBar.setProgress,
        		self.binningProgressBar.setTotalSteps, self.binningProgressBar.setProgress)
        	datafile=open(str(self.dataFileEdit.text()),'w')
        	for i in self.chrhist.keys():
        		print >> datafile, "#%s" % (str(self.geo.PTE[i]))
        		for j in range(len(self.chrhist[i][0])):
        			print >> datafile, "%e   %e" % (self.chrhist[i][0][j],self.chrhist[i][1][j])
        		print >> datafile, '\n'
        	datafile.close()
        	if self.callGnuplotCheckBox.isChecked():
        		self.plotHistogram()
        	QDialog.accept(self)
        

    def plotHistogram(self):
        	chrplot=[]
        	for i in self.chrhist.keys():
        		chrplot.append(self.Gnuplot.Data(self.chrhist[i][0],self.chrhist[i][1],title=self.geo.PTE[i],inline=1))
        	self.plot.plot(*chrplot)
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
        		self.plot.plot(*chrplot)
        		self.plot("set out")
        

    def setGeometry(self,a0):
        	self.geo=a0
        	self.setFinishEnabled(self.WizardPage_2,True)
        

    def init(self):
        	"""Initialize Charge histogram wizard
        	check for Gnuplot package and init that if available."""
        	try:
        		self.Gnuplot = __import__("Gnuplot")
        	except:
        		self.callGnuplotCheckBox.setEnabled(False)
        	else:
        		self.callGnuplotCheckBox.setEnabled(True)
        		self.plot=self.Gnuplot.Gnuplot(debug=0,persist=1)
        		self.plot("set data style boxes")
        		self.plot("set style fill solid")
        		self.plot("set title 'element atom charges histograms'")
        		self.plot("set xlabel 'atomic charge [e^{-}]'")
        		self.plot("set ylabel '# of atoms'")
        

    def __tr(self,s,c = None):
        return qApp.translate("elementChargesHistograms",s,c)
