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


void atomChargeConstraintsWizard::printConstraints()
{
	# parse the atom serial numbers, get the DFTB+ charge constraints string, put them into the output textEdit and switch pages
	serials=str(self.atomSerialsTextEdit.text())
	# user interface counts from 1, internal indices count from 0
	try:
		indices=[ int(j)-1 for j in serials.split() ]
	except ValueError:
		QMessageBox.warning(self,"Invalid atom numbers","<H1>Error</H1><p>The list of atom serials could not be understood. Atom serials must be given as a white-space separated list of integers. Please check your input.</p>")
	# convert prefactor to float
	prefactor=float(self.constraintsPrefactorLineEdit.text())
	try:
		self.constraintsDisplayTextEdit.setText(self.geo.getHSDChargeConstraintString(indices, prefactor))
	except ValueError:
		QMessageBox.warning(self,"Invalid atom number","<H1>Error</H1><p>One or more atom serials are not present in the geometry. Please check your input.</p>")
	except:
		raise
	# enable finish button
	self.setFinishEnabled(self.WizardPage,True)
}


void atomChargeConstraintsWizard::init()
{
	self.constraintsPrefactorLineEdit.setText("1.00")
}


void atomChargeConstraintsWizard::pageFunctionDispatch()
{

}


void atomChargeConstraintsWizardBase::setGeometry( geo )
{
	# store a link to the geometry object to work on
	self.geo=a0
}
