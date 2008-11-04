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


void periodicExpand::setVectors( Lattice )
{
	Lattice=a0*5.291772e-01
	self.aVectorLabel.setText("&lt;%7.4f %7.4f %7.4f&gt;" % (Lattice[0][0],Lattice[0][1],Lattice[0][2]))
	self.bVectorLabel.setText("&lt;%7.4f %7.4f %7.4f&gt;" % (Lattice[1][0],Lattice[1][1],Lattice[1][2]))
	self.cVectorLabel.setText("&lt;%7.4f %7.4f %7.4f&gt;" % (Lattice[2][0],Lattice[2][1],Lattice[2][2]))
}


void periodicExpand::getExpandSpecifier()
{
	return (self.aSpinBox.value(),self.bSpinBox.value(),self.cSpinBox.value())
}
