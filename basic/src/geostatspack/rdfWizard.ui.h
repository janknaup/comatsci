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
	self.rdfdata=self.geo.rdf(self.stepwidth, 
		self.postProcessProgressBar.setTotalSteps, self.postProcessProgressBar.setProgress,
		self.binningProgressBar.setTotalSteps, self.binningProgressBar.setProgress)
	datafile=open(str(self.dataFileEdit.text()),'w')
	print >> datafile, "# r[a.u.]  rdf"
	for j in range(len(self.rdfdata[0])):
		print >> datafile, "%e   %e" % (self.rdfdata[0][j],self.rdfdata[1][j])
	datafile.close()
	if self.callGnuplotCheckBox.isChecked():
		self.plotHistogram()
	QDialog.accept(self)
}


void rdfWizard::plotHistogram()
{
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
		self.plot("set title 'total radial distribution function'")
		self.plot("set xlabel 'r[Bohr]'")
		self.plot("set ylabel 'rdf'")
	#initialize stepwidth input slider and LCD display
	self.stepWidthSlider.setValue(20)
	self.stepWidthLCDNumber.display(0.20)
	self.stepwidth=0.20
}


void rdfWizard::stepWidthChanged(int)
{
	self.stepwidth=float(a0)/100.0
	self.stepWidthLCDNumber.display(self.stepwidth)
}
