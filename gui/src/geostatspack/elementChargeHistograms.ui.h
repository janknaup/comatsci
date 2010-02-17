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


void elementChargesHistograms::setSavePlotState(bool)
{
	self.savePlotCheckBox.setEnabled(a0)
	self.plotFileEdit.setEnabled(a0)
	self.plotFilenameTextLabel.setEnabled(a0)
	self.fileFormatButtonGroup.setEnabled(a0)
}



void elementChargesHistograms::accept()
{
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
}


void elementChargesHistograms::plotHistogram()
{
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
}


void elementChargesHistograms::setGeometry(geo)
{
	self.geo=a0
	self.setFinishEnabled(self.WizardPage_2,True)
}


void elementChargesHistograms::init()
{
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
}
