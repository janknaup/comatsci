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


void layerLinkBCTC::setGeometry(geo)
{
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
		}





void layerLinkBCTC::my_link()
{
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
}


void layerLinkBCTC::setBCTCTable(BCTcoefficients)
{
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
}


void layerLinkBCTC::BCTCCoefficientsFromTable()
{
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
}


void layerLinkBCTC::BCTCCoefficientsFromFile()
{
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
}


void layerLinkBCTC::BCTCCoefficientsFromCalc()
{
	# calculate BCTC coefficients from geometry and apply to the coefficients table
	BCTCoefficients=self.geo.getElementElementChargeTransfers()[0]
	self.setBCTCTable(BCTCoefficients)
}


void layerLinkBCTC::linkAction()
{
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
}
