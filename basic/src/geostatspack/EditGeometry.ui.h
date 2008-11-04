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

void Edit_Atoms::addLayer()
{
	lname=str(self.layerNameEdit.text())
	li=self.geo.addlayer(lname)
	newrow=self.layerTable.numRows()
	self.layerTable.setNumRows(newrow+1)
	self.layerTable.setItem(newrow,0,QTableItem(self.layerTable,QTableItem.Never,str(li)))
	self.layerTable.setItem(newrow,1,QTableItem(self.layerTable,QTableItem.Never,lname))
	layerlabelstring="%i: %s" % (li,lname)
	self.layersqsl.append(layerlabelstring)
	self.updateLayersInAtomTable()
	self.layerNameEdit.setText("")
		}

void Edit_Atoms::my_ok()
{
	self.my_apply()
	self.close()
}


void Edit_Atoms::setGeometry( geo )
{
	#keep a the original geometry untouched
	self.originalgeo=a0
	self.geo=copy.deepcopy(a0)
	#initialize layer table and Layer Spinboxes
	layercount=len(self.geo.LayerDict)
	self.layerTable.setNumRows(layercount)
	self.layersqsl=QStringList()
	for i in range(layercount):
		li=self.geo.LayerDict.keys()[i]
		lname=self.geo.LayerDict[li].Name
		layerlabelstring="%i: %s" % (li,lname)
		self.layerTable.setItem(i,0,QTableItem(self.layerTable,QTableItem.Never,str(li)))
		self.layerTable.setItem(i,1,QTableItem(self.layerTable,QTableItem.Never,str(lname)))
		self.layersqsl.append(layerlabelstring)
	self.layerTable.adjustColumn(0)
	self.layerTable.setColumnStretchable(1,True)
	#initialize atom table
	self.PTEqsl=QStringList()
	for i in self.geo.PTE:
		self.PTEqsl.append(i)
	self.atomTable.setNumRows(self.geo.Atomcount)
	for i in range(self.geo.Atomcount):
		self.atomTable.setItem(i,0,QComboTableItem(self.atomTable,self.PTEqsl,False))
		self.atomTable.item(i,0).setCurrentItem(self.geo.AtomTypes[i])
		self.atomTable.setItem(i,1,QTableItem(self.atomTable,QTableItem.WhenCurrent,str(self.geo.Geometry[i][0])))
		self.atomTable.setItem(i,2,QTableItem(self.atomTable,QTableItem.WhenCurrent,str(self.geo.Geometry[i][1])))
		self.atomTable.setItem(i,3,QTableItem(self.atomTable,QTableItem.WhenCurrent,str(self.geo.Geometry[i][2])))
		self.atomTable.setItem(i,4,QComboTableItem(self.atomTable,self.layersqsl,False))
		self.atomTable.item(i,4).setCurrentItem(self.geo.AtomLayers[i])
		self.atomTable.setItem(i,5,QTableItem(self.atomTable,QTableItem.WhenCurrent,str(self.geo.AtomCharges[i])))
		self.atomTable.setItem(i,6,QTableItem(self.atomTable,QTableItem.WhenCurrent,str(self.geo.AtomSubTypes[i])))
	for i in range(6):
		self.atomTable.adjustColumn(i)
	self.atomTable.setColumnStretchable(6,True)
	}
	


void Edit_Atoms::deleteAtom()
{
	curatom=self.atomTable.currentRow()
	self.atomTable.removeRow(curatom)
}


void Edit_Atoms::revert()
{
	self.setGeometry(self.originalgeo)
}


void Edit_Atoms::my_apply()
{
	#construct a  completely new geometry object for output
	newgeo=self.geo.__class__()
	#first set the new geometry's Layer Dictionary     '
	for i in range(self.layerTable.numRows()):
		li=int(str(self.layerTable.text(i,0)))
		lname=str(self.layerTable.text(i,1))
		newgeo.LayerDict[li]=Geometry.GeoLayer(lname)
	#just reuse mode Lattice and origin
	newgeo.Origin=self.geo.Origin
	newgeo.Lattice=self.geo.Lattice
	newgeo.Mode=self.geo.Mode
	#now populate the new Geometry with atoms
	for i in range(self.atomTable.numRows()):
		type=self.atomTable.item(i,0).currentItem()
		x=float(str(self.atomTable.text(i,1)))
		y=float(str(self.atomTable.text(i,2)))
		z=float(str(self.atomTable.text(i,3)))
		li=self.atomTable.item(i,4).currentItem()
		chr=float(str(self.atomTable.text(i,5)))
		st=str(self.atomTable.text(i,6))
		newgeo.addatom(type,(x,y,z),li,chr,st)
	#This cunningly overwrites the original geometry object
	self.setGeometry(newgeo)
}


void Edit_Atoms::addAtom()
{
	curatom=self.atomTable.currentRow()
	type=self.atomTable.item(curatom,0).currentItem()
	li=self.atomTable.item(curatom,4).currentItem()
	chr=self.atomTable.text(curatom,5)
	st=self.atomTable.text(curatom,6)
	self.atomTable.insertRows(curatom)
	self.atomTable.setItem(curatom,0,QComboTableItem(self.atomTable,self.PTEqsl,False))
	self.atomTable.item(curatom,0).setCurrentItem(type)
	self.atomTable.setText(curatom,1,"")
	self.atomTable.setText(curatom,2,"")
	self.atomTable.setText(curatom,3,"")
	self.atomTable.setItem(curatom,4,QComboTableItem(self.atomTable,self.layersqsl,False))
	self.atomTable.item(curatom,4).setCurrentItem(li)
	self.atomTable.setText(curatom,5,"0.0")
	self.atomTable.setText(curatom,6,st)
}


void Edit_Atoms::updateLayersInAtomTable()
{
	for i in range(self.atomTable.numRows()):
		olditem=self.atomTable.item(i,4).currentItem()
		self.atomTable.setItem(i,4,QComboTableItem(self.atomTable,self.layersqsl,False))
		self.atomTable.item(i,4).setCurrentItem(olditem)
}


void Edit_Atoms::deleteLayer()
{
	curlayer=self.layerTable.currentRow()
	# Do not touch the default layer
	if curlayer!=0:
		self.layerTable.removeRow(curlayer)
		#renumber the layers
		for i in range(curlayer-1,self.layerTable.numRows()):
			self.layerTable.setText(i,0,str(i))
		#put atoms from deleted layer into default layer
		#and make atoms with layer indices higher than deleted layer consistent again
		for i in range(self.atomTable.numRows()):
			itmlr=self.atomTable.item(i,4).currentItem()
			if itmlr==curlayer:
				self.atomTable.item(i,4).setCurrentItem(0)
			elif itmlr>curlayer:
				self.atomTable.item(i,4).setCurrentItem(itmlr-1)
		#make layers QStringList consistent with table and update comboBoxes for all atoms
		del self.layersqsl[curlayer]
		self.updateLayersInAtomTable()
}
