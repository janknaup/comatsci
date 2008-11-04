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


void layerLinkBase::setGeometry(geo)
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
		}





void layerLinkBase::my_link()
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
