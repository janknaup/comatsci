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
	serials=self.atomSerialsTextEdit.getText()
	# user interface counts from 1, internal indices count from 0
	indices=[ int(j)-1 for j in serials.split() ]
	# convert prefactor to float
	prefactor=float(self.atomSerialsTextEdit.getText())
	self.constraintsDisplayTextEdit.setText(self.geo.getHSDChargeConstraintString(indices, prefactor))
	# enable finish button
	self.setFinishEnabled(True)
}


void atomChargeConstraintsWizard::init()
{
	self.constraintsPrefactorLineEdit.setText("1.00")
}


void atomChargeConstraintsWizard::pageFunctionDispatch(QString)
{
# call page-specific functions to initialize values
if a0=="ChargeConstraints":
	self.printConstraints()
}
