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


void Form1::fileOpen()
{
	fn=str(QFileDialog.getOpenFileName(".","Known formats (*.gen *.xyz *.fmg)", self))
	if (len(fn)!=0):
		self.fileread(fn)
}



void Form1::fileSave()
{
	if self.filename.endswith(".gen"):
		self.geo.writegen(self.filename)
	elif self.filename.endswith(".xyz"):
		self.geo.writexyz(self.filename)
	elif self.filename.endswith(".fdf"):
		self.geo.writefdf(self.filename)
	elif self.filename.endswith(".pdb"):
		self.geo.writepdb(self.filename)
	elif self.filename.endswith(".fmg"):
		self.geo.writefmg(self.filename)
}


void Form1::fileSaveAs()
{
	fd=QFileDialog(".","Generic (*.gen)",self)
	fd.setMode(QFileDialog.AnyFile)
	fd.addFilter("XYZ (*.xyz)")
	fd.addFilter("SIESTA (*.fdf)")
	fd.addFilter("PDB (*.pdb)")
	fd.addFilter("FMG (*.fmg)")
	fd.setSelectedFilter("gen")
	if (fd.exec_loop()==QDialog.Accepted):
		filename=str(fd.selectedFile())
		if filename.endswith(".gen"):
			self.geo.writegen(filename)
		elif filename.endswith(".xyz"):
			self.geo.writexyz(filename)
		elif filename.endswith(".fdf"):
			self.geo.writefdf(filename)
		elif filename.endswith(".pdb"):
			self.geo.writepdb(filename)
		elif filename.endswith(".fmg"):
			self.geo.writefmg(filename)
}


void Form1::DisplayCoordinations()
{
	self.textBrowser1.append(self.geo.rt_coordinations())
	self.statisticsCoordinationsAction.setDisabled(1)
}



void Form1::helpContents()
{
	gsh=geostatshelp.GeostatsHelp()
	gsh.show()
}


void Form1::helpAbout()
{
	QMessageBox.about(self,"About geostats","""<H1>geostats %s</H1><p>Calculate geometry statistics and transform file formats.</p>
			  <p>Licensed under the Non-Profit Open Software License version 3.0</p>
			  <p>see file LICENSE for details.</p>
			  <p>Written by Jan M. Knaup <b><a href="mailto:Jan.Knaup@bccms.uni-bremen.de">Jan.Knaup@bccms.uni-bremen.de</a></b></p>"""%(comatsci.constants.VERSION))
}



void MainWindow::SaveStatistics()
{
	statsfile=open("geostats.htm","w")
	statsfile.write(str(self.textBrowser1.text()))
}


void MainWindow::fileread(fnam)
{
	#save filename for later, initialize geometry object and read file
	self.filename=a0
	self.geo=comatsci.Geometry.FullFeaturedGeometry()
	self.geo.readfile(a0)
	#prepare the statistics text display
	statsheader="<H1>Statistics on %s</H1>" % (a0.rsplit("/",1)[1])
	self.textBrowser1.setText(statsheader)
	#enable general operations
	self.fileSaveAsAction.setEnabled(True)
	self.fileSaveAction.setEnabled(True)
	self.editAtomsButton.setEnabled(True)
	self.editedit_atomsAction.setEnabled(True)
	#enable linkatoms generation if more than one layer exists
	if len(self.geo.LayerDict)>1 and ("embed" in self.geo.getFeatures()):
		self.editgenerate_linkatomsAction.setEnabled(True)
		self.editHCS_link_layersAction.setEnabled(True)
		self.editBCTC_link_LayersAction.setEnabled(True)
##	self.VMDButton.setEnabled(True)
##	self.fileView_in_VMDAction.setEnabled(True)
	#only enable periodic expansion for supercell geometries
	if self.geo.Mode=="S":
		self.editperiodic_expandAction.setEnabled(True)
		self.periodicExpandButton.setEnabled(True)
	#enable statistics actions after file read
	self.statisticsbond_listAction.setEnabled(True)
	self.statisticsChargesAction.setEnabled(True)
	self.statisticsCoordinationsAction.setEnabled(True)
	self.statisticsgraphsAction.setEnabled(True)
	self.statisticsgraphselement_charge_histogramsAction.setEnabled(True)
	self.statisticsgraphsradial_distribution_functionAction.setEnabled(True)
	self.statisticsSave_StatisticsAction.setEnabled(True)
	self.statisticsget_charge_constraintsAction.setEnabled(True)
	self.statisticssave_BCT_coefficientsAction.setEnabled(True)
}


void MainWindow::display_bondlist()
{
	self.textBrowser1.append(self.geo.rt_bondlist())
	self.statisticsbond_listAction.setDisabled(1)
}


void MainWindow::editAtoms()
{
	GeoEditor=EditGeometry.Edit_Atoms()
	GeoEditor.setGeometry(self.geo)
	GeoEditor.exec_loop()
	self.geo=GeoEditor.originalgeo
	if self.geo.layerbyname("PCHR")!=None:
		self.editgenerate_linkatomsAction.setEnabled(True)
}


void MainWindow::periodicExpand()
{
	pedialog=periodicExpand.periodicExpand()
	pedialog.setVectors(self.geo.Lattice)
	if pedialog.exec_loop()==QDialog.Accepted:
		self.geo.periodicexpand(pedialog.getExpandSpecifier())
}


void MainWindow::DisplayCharges()
{
	self.textBrowser1.append(self.geo.rt_charges())
	self.statisticsChargesAction.setDisabled(1)
}


void MainWindow::generate_linkatoms()
{
	linkDialog=linkLayersSLA.layerLinkSLA()
	linkDialog.setGeometry(self.geo)
	if linkDialog.exec_loop()==QDialog.Accepted:
		self.geo=linkDialog.embeddedGeometry
}

void MainWindow::generate_hcsLinkAtoms()
{
	linkDialog=linkLayersHCS.layerLinkHCS()
	linkDialog.setGeometry(self.geo)
	if linkDialog.exec_loop()==QDialog.Accepted:
		self.geo=linkDialog.embeddedGeometry
}

void MainWindow::elementChargesHistograms()
{
	elcharhist=elementChargeHistograms.elementChargesHistograms()
	elcharhist.setGeometry(self.geo)
	elcharhist.exec_loop()
}


void MainWindow::rdf()
{
	rdfwiz=rdfWizard.rdfWizard()
	rdfwiz.setGeometry(self.geo)
	rdfwiz.exec_loop()
}


void MainWindow::generate_BCTCLinkAtoms()
{
	linkDialog=linkLayersBCTC.layerLinkBCTC()
	linkDialog.setGeometry(self.geo)
	if linkDialog.exec_loop()==QDialog.Accepted:
		self.geo=linkDialog.embeddedGeometry
}


void MainWindow::save_BCTC_coefficients()
{
	# pop-up a file selection dialog
	fileDialog=QFileDialog(".","whitespace separated data (*.bctc)",self)
	fileDialog.setMode(QFileDialog.AnyFile)
	
	# if the user didnt abort, write to the given file
	if (fileDialog.exec_loop()==QDialog.Accepted):
		filename=str(fileDialog.selectedFile())
		# check if filename ends with .dat, append proper extension if not
		if not filename[-4:-1] in [".bctc",".BCTC"]:
			filename+=".bctc"
		BCTCFile=open(filename,'w')
		# write a comment line describing the columns
		print >> BCTCFile, "#Z1   Z2   charge transfer"
		# get BCTC coeffcients dictionary
		coefficients=self.geo.getElementElementChargeTransfers()[0]
		# iterate through keys. The key is a 2-tuple of nuclear charges
		for i in coefficients.keys():
			print >> BCTCFile, "%2d   %2d   %f" % (i[0],i[1],coefficients[i])
		# finished
		BCTCFile.close()
}


void MainWindow::display_charge_constraints()
{
	# we have a wizard for this, just for convenience
	constraintsDialog=chargeConstraintsWizard.atomChargeConstraintsWizard()
	# pass the geomerty object to the wizard
	constraintsDialog.setGeometry(self.geo)
	# The wizard does everything by itself, so just run it
	constraintsDialog.exec_loop()
}
