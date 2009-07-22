/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you want to add, delete, or rename functions or slots, use
** Qt Designer to update this file, preserving your code.
**
** You should not define a constructor or destructor in this file.
** Instead, write your code in functions called init() and destroy().
** These will automatically be called by the form's constructor and
** destructor.
*****************************************************************************/


void rdfWizard::setSavePlotState( bool )
{
	self.savePlotCheckBox.setEnabled(a0)
	self.plotFileEdit.setEnabled(a0)
	self.plotFilenameTextLabel.setEnabled(a0)
	self.fileFormatButtonGroup.setEnabled(a0)
}


void rdfWizard::accept()
{
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
}


void rdfWizard::plotHistogram()
{
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
}


void rdfWizard::setGeometry( geo )
{
	self.geo=a0
	self.setFinishEnabled(self.WizardPage_2,True)
}


void rdfWizard::init()
{
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
}


void rdfWizard::binsChanged(int)
{
	self.bins=int(a0)
	self.binsLCDNumber.display(self.bins)
}
