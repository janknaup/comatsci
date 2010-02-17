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


void bondlengthsWizard::setSavePlotState( bool )
{
	self.savePlotCheckBox.setEnabled(a0)
	self.plotFileEdit.setEnabled(a0)
	self.plotFilenameTextLabel.setEnabled(a0)
	self.fileFormatButtonGroup.setEnabled(a0)
}


void bondlengthsWizard::accept()
{
	self.blHist=self.geo.getBondLengthHistograms(self.bins)
	keys=self.blHist.keys()
	keys.remove((-1,-1))
	datafile=open(str(self.dataFileEdit.text()),'w')
	print >> datafile, "# r[a.u.] \ttotal", 
	for i in keys:
		print >> datafile, "\t%s" %(self.geo.PTE[i[0]]+"-"+self.geo.PTE[i[1]]),
	print >> datafile, ""
	for j in range(len(self.blHist[(-1,-1)][0])):
		print >> datafile, "%e\t%e" % (self.blHist[(-1,-1)][0][j],self.blHist[(-1,-1)][1][j]),
		for i in keys:
			print >> datafile, "\t%e" % self.blHist[i][1][j],
		print >> datafile, ""
	datafile.close()
	if self.callGnuplotCheckBox.isChecked():
		self.plotHistogram()
	QDialog.accept(self)
}


void bondlengthsWizard::plotHistogram()
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


void bondlengthsWizard::setGeometry( geo )
{
	self.geo=a0
	self.setFinishEnabled(self.WizardPage_2,True)
}


void bondlengthsWizard::init()
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
		self.plot("set title 'total bond length histogram'")
		self.plot("set xlabel 'r[Bohr]'")
		self.plot("set ylabel 'count'")
	#initialize bins input slider and LCD display
	self.binsSlider.setValue(40)
	self.binsLCDNumber.display(40)
	self.bins=40
}


void bondlengthsWizard::binsChanged(int)
{
	self.bins=int(a0)
	self.binsLCDNumber.display(self.bins)
}
