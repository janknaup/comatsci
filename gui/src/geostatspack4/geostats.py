#!/usr/bin/python
##############################################################################
# geostats.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-20012 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################
from __future__ import print_function

import sys,codecs
import comatsci
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import copy

#from optparse import OptionParser

from PyQt4 import QtGui,QtCore

import geostatspack4

class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.).
    Based upon Matplotlib Qt4 example
         Copyright (C) 2005 Florent Rougon
                       2006 Darren Dale"""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes to remain
        self.axes.hold(True)

        self.compute_initial_figure()

        #
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass



class MainWindow(QtGui.QMainWindow):
    
    saveModes={"xyz" : ("xyz",comatsci.geometry.Geometry.writexyz),
               "pdb" : ("pdb",comatsci.geometry.Geometry.writepdb),
               "fdf" : ("fdf",comatsci.geometry.Geometry.writefdf),
               "gen" : ("gen",comatsci.geometry.Geometry.writegen),
               "fmg" : ("fmg",comatsci.geometry.Geometry.writefmg),
               "xyzq": ("xyzq",comatsci.geometry.Geometry.writexyzq),
               "tm"  : ("coord",comatsci.geometry.Geometry.writeTurboMole),
               "aims": ("geometry.in",comatsci.geometry.Geometry.writeAIMS),
               "cdh" : ("cdh",comatsci.geometry.Geometry.writeCDH)
               }
    
    fixedFileNameTypes=("aims","tm")
    filterFileTypes={"xmol xyz (*.xyz)"                     :   "xyz" ,
                     "protein database format (*.pdb)"      :   "pdb" ,
                     "SIESTA flexible data format (*.fdf)"  :   "fdf" ,
                     "DFTB generic (*.gen)"                 :   "gen", 
                     "flexible molecular geometry (*.fmg)"  :   "fmg" ,
                     "DFTB point charge format (*.xyzq)"    :   "xyzq",
                     "TurboMole geometry format (coord)"    :   "tm"  ,
                     "FHI aims input (geometry.in)"         :   "aims",
                     "chemical data hierarchy (*.cdh)"      :   "cdh" 
                     }
    
    def __init__(self,parent=None):
        # setup widget and user interface
        QtGui.QWidget.__init__(self, parent)
        self.ui = geostatspack4.geostatsMain.Ui_MainWindow()
        self.ui.setupUi(self)
        
        # connect all the menu signals
        QtCore.QObject.connect(self.ui.actionHelp, QtCore.SIGNAL("triggered()"),self.showHelp)
        QtCore.QObject.connect(self.ui.actionAbout, QtCore.SIGNAL("triggered()"),self.showAbout)
        QtCore.QObject.connect(self.ui.actionPrint_Output, QtCore.SIGNAL("triggered()"),self.handlePrint)
        QtCore.QObject.connect(self.ui.actionSave_output, QtCore.SIGNAL("triggered()"),self.handleSave_Output)
        QtCore.QObject.connect(self.ui.actionOpen, QtCore.SIGNAL("triggered()"),self.handleOpen)
        QtCore.QObject.connect(self.ui.actionCoordinations, QtCore.SIGNAL("triggered()"),self.handleCoordinations)
        QtCore.QObject.connect(self.ui.actionCharges, QtCore.SIGNAL("triggered()"),self.handleCharges)
        QtCore.QObject.connect(self.ui.actionSave_geometry, QtCore.SIGNAL("triggered()"),self.handleSave_geometry)
        QtCore.QObject.connect(self.ui.actionSave_geometry_as, QtCore.SIGNAL("triggered()"),self.handleSave_geometry_as)
        QtCore.QObject.connect(self.ui.actionPeriodic_Expand, QtCore.SIGNAL("triggered()"),self.handlePeriodicExpand)
        QtCore.QObject.connect(self.ui.actionRadial_Distribution_Functions, QtCore.SIGNAL("triggered()"),self.handleRadialDistributionFunctions)
        QtCore.QObject.connect(self.ui.actionElement_Charge_Histograms, QtCore.SIGNAL("triggered()"),self.handleElementChargeHistograms)
        QtCore.QObject.connect(self.ui.actionBond_Length_Histogram, QtCore.SIGNAL("triggered()"),self.handleBondLengthHistograms)
        QtCore.QObject.connect(self.ui.actionAngle_Distribution_Histogram, QtCore.SIGNAL("triggered()"),self.handleBondAngleHistograms)
        QtCore.QObject.connect(self.ui.actionDFTB_Charge_Constraints, QtCore.SIGNAL("triggered()"),self.handleChargeConstraint)
        QtCore.QObject.connect(self.ui.actionSimple_link_Layers, QtCore.SIGNAL("triggered()"),self.handleSimpleLink)
        QtCore.QObject.connect(self.ui.actionBCTC_link_Layers, QtCore.SIGNAL("triggered()"),self.handleBCTCLink)
        QtCore.QObject.connect(self.ui.actionWrite_BCTC_Coefficients, QtCore.SIGNAL("triggered()"),self.handleSaveBCTC)
        QtCore.QObject.connect(self.ui.actionEdit_Atoms, QtCore.SIGNAL("triggered()"),self.handleEditAtoms)
        
        
        
    #*****************************************************************************
    #  Actual application logic  
    #*****************************************************************************
    def showHelp(self):
        """Display help dialog
        """
        dlghelp=GeostatsHelp()
        helpTextFile = QtCore.QFile(":/html/html/geostatsHelp.htm")
        if helpTextFile.open(QtCore.QIODevice.ReadOnly):
            helpString=QtCore.QString(helpTextFile.readAll())
            dlghelp.ui.textBrowser.setHtml(unicode(helpString).encode("ascii",'xmlcharrefreplace'))
            dlghelp.ui.textBrowser.scrollToAnchor("head")
        dlghelp.exec_()

        
    def showAbout(self):
        about=QtGui.QMessageBox(parent=self)
        about.setText("About geostats.py")
        infoTextFile = QtCore.QFile(":/html/html/about.htm")
        if infoTextFile.open(QtCore.QIODevice.ReadOnly):
            infoString=QtCore.QString(infoTextFile.readAll())
        else:
            infoString="part of the comatsci computational materials science toolkit\nCopyright (c) 2004-2012 Jan M. Knaup\nError reading addition information"
        about.setInformativeText(infoString)
        infoTextFile.close()
        about.setIconPixmap(QtGui.QPixmap(":/images/images/GSLogo.svg"))
        about.show()
    
    
    def handlePrint(self):
        dialog = QtGui.QPrintDialog()
        if dialog.exec_() == QtGui.QDialog.Accepted:
            self.ui.textBrowser.document().print_(dialog.printer())
    
            
    def handleSaveBCTC(self):
        bctcFileName=str(QtGui.QFileDialog.getSaveFileName(parent=self, caption=QtCore.QString("Open BCTC Coefficients"),
                                                      filter=QtCore.QString("BCTC coefficient files (*.bctc);;all supported files (*.*)")))
        if bctcFileName!="":
            if not bctcFileName.endswith(".bctc"):
                bctcFileName+=".bctc"
            tempCTC=self.geo.getElementElementChargeTransfers()[0]
            bctcFile=open(bctcFileName,"w")
            print("#Elements  Bond_Charge_Transfer_Coefficient",file=bctcFile)
            for coefficient in tempCTC.keys():
                print("{0:d}\t{1:d}\t{2:10.6f}".format(coefficient[0],coefficient[1],tempCTC[coefficient]),file=bctcFile)
            bctcFile.close()
        
    
    
    def handleSave_Output(self):
        (savefile, savefilter)=QtGui.QFileDialog.getSaveFileNameAndFilter(parent=self, caption=QtCore.QString("Save Output"),
                                                                          filter=QtCore.QString("plain text (*.txt);;html (*.htm *.html)"))#;;PDF (*.pdf)"))
        if savefilter=="plain text (*.txt)":
            if savefile.toLower()[-4:]!=".txt":
                savefile+=".txt"
            outfile=codecs.open(savefile,"w",encoding="ascii")
            print(unicode(self.ui.textBrowser.document().toPlainText()).encode("ascii",'replace'),file=outfile)
            outfile.close()
        elif savefilter=="html (*.htm *.html)":
            if savefile.toLower()[-4:]!=".htm" and savefile.toLower()[-5:]!=".html":
                savefile+=".html"
            outfile=open(savefile,"w")
            print(unicode(self.ui.textBrowser.document().toHtml(encoding='utf-8')).encode("ascii",'xmlcharrefreplace'),file=outfile)
            outfile.close()
#        elif savefilter=="PDF (*.pdf)":
#            if savefile.toLower()[-4:]!=".pdf":
#                savefile+=".pdf"
#            pass
        else:
            pass
    
    def handleSimpleLink(self):
        embedDialog=SimpleLinkDialog(self.geo, self)
        if embedDialog.exec_():
            self.geo = embedDialog.linkedGeo
            self.resetForGeoChanged("SLAlinked.fmg")
            
    
    def handleBCTCLink(self):
        embedDialog=BCTCLinkDialog(self.geo, self)
        if embedDialog.exec_():
            self.geo = embedDialog.linkedGeo
            self.resetForGeoChanged("BCTClinked.fmg")
        
    
    def handleOpen(self):
        geoFileName=QtGui.QFileDialog.getOpenFileName(parent=self, caption=QtCore.QString("Open Geometry"),
                                                      filter=QtCore.QString("all supported files (*.gen *.xyz *.cdh *.fmg geometry.in);;DFTB generic (*.gen);;xmol xyz (*.xyz);;chemical data hierarchy (*.cdh);;FHI aims input (geometry.in);;flexible molecular geometry (*.fmg)"))
        if not geoFileName.isEmpty():
            self.readGeometry(str(geoFileName))        
        else:
            pass
        
    def handleSave_geometry(self):
        filemode=self.geoFileName.lower()[self.geoFileName.lower().rfind(".")+1:]
        if filemode!=-1:
            self.saveModes[filemode][1](self.geo,self.geoFileName)
        else:
            errorDialog=QtGui.QErrorMessage(parent=self)
            errorDialog.showMessage("Error determining file type from file name '{0:s}'".format(self.geoFileName))
    
    
    def handleSave_geometry_as(self):
        filenamePrefix=self.geoFileName.lower()[0:self.geoFileName.lower().rfind(".")]
        print(filenamePrefix)
        (savefile, savefilter)=QtGui.QFileDialog.getSaveFileNameAndFilter(parent=self, caption=QtCore.QString("Save Geometry As"),
                                                                          filter=QtCore.QString(";;".join(self.filterFileTypes.keys())))
        savefile=str(savefile)
        savefilter=str(savefilter)
        print(savefile,savefilter)
        if self.filterFileTypes.has_key(savefilter):
            filenamesuffix=self.saveModes[self.filterFileTypes[savefilter]][0]
            if self.filterFileTypes[savefilter] in self.fixedFileNameTypes:
                savefile=filenamesuffix
                errorDialog=QtGui.QErrorMessage(parent=self)
                errorDialog.setModal(True)
                errorDialog.showMessage("Attention, filename changed to '{0:s}'".format(savefile))
            else:
                if savefile.find(filenamesuffix)==-1:
                    savefile=savefile+"."+filenamesuffix
            try:
                self.saveModes[self.filterFileTypes[savefilter]][1](self.geo,savefile)
            except comatsci.geometry.GeometryError as error:
                errorDialog=QtGui.QErrorMessage(parent=self)
                errorDialog.setModal(True)
                errorDialog.showMessage(error.args[0])
            
        
    

    def resetForGeoChanged(self,geoFileName):
        self.ui.textBrowser.document().clear()
        self.ui.actionSave_geometry.setEnabled(True)
        self.ui.actionSave_geometry_as.setEnabled(True)
        self.ui.menuSummaries.setEnabled(True)
        self.ui.actionCharges.setEnabled(True)
        self.ui.actionCoordinations.setEnabled(True)
        self.ui.menuGraphs.setEnabled(True)
        self.ui.actionRadial_Distribution_Functions.setEnabled(True)
        self.ui.actionWrite_BCTC_Coefficients.setEnabled(True)
        self.ui.actionSimple_link_Layers.setEnabled(True)
        self.ui.actionBCTC_link_Layers.setEnabled(True)
        self.ui.actionDFTB_Charge_Constraints.setEnabled(True)
        self.ui.actionEdit_Atoms.setEnabled(True)
        if self.geo.Mode=="S":
            self.ui.actionPeriodic_Expand.setEnabled(True)
        self.ui.textBrowser.insertHtml("<H1>Geometry Statistics</H1>\n<p>statistics for file &lsquo;{0:s}&rsquo;</p>\n<br>\n".format(geoFileName))
    
    
    
    def handleChargeConstraint(self):
        constraintDialog=chargeConstraintsDialog(self.geo, self)
        constraintDialog.exec_()
        


    def readGeometry(self,geoFileName):
        self.geoFileName=str(geoFileName)
        self.geo=comatsci.geometry.qmmmGeometry()
        self.geo.readfile(str(geoFileName))
        self.resetForGeoChanged(geoFileName)
        

    
    def handleCoordinations(self):
        coordinationsSummary=self.geo.rt_coordinations()
        self.ui.textBrowser.insertHtml(coordinationsSummary)
        self.ui.actionCoordinations.setEnabled(False)
    
    
    def handleCharges(self):
        chargesSummary=self.geo.rt_charges()
        self.ui.textBrowser.insertHtml(chargesSummary)
        self.ui.actionCharges.setEnabled(False)
        
    
    def handlePeriodicExpand(self):
        expandDialog=PeriodicExpand(parent=self)
        result=expandDialog.exec_()
        if result:
            aRepeat=expandDialog.ui.aSpinBox.value()
            bRepeat=expandDialog.ui.bSpinBox.value()
            cRepeat=expandDialog.ui.cSpinBox.value()
            self.geo.periodicexpand((aRepeat,bRepeat,cRepeat))
            self.resetForGeoChanged(self.geoFileName)
            
            
    def handleEditAtoms(self):
        editDialog=EditAtomsDialog(self.geo,parent=self)
        result=editDialog.exec_()
        if result:
            self.geo=editDialog.geo
            self.resetForGeoChanged(self.geoFileName)
            
    
    def handleRadialDistributionFunctions(self):
        plotDialog=RDFDialog(geometry=self.geo,parent=self)
        result=plotDialog.exec_()
        if result and plotDialog.ui.addOutputCheckBox.isChecked():
            plotDialog.ui.plotCanvas.print_svg("rdfplot.svg")
            self.ui.textBrowser.insertHtml("<H2>Radial Pair Distribution Function<H2>\n")
            self.ui.textBrowser.insertHtml('<p><img src="rdfplot.svg" /></p><br />\n')
    
    
    def handleElementChargeHistograms(self):
        plotDialog=ElementChargeHistogramDialog(geometry=self.geo,parent=self)
        result=plotDialog.exec_()
        if result and plotDialog.ui.addOutputCheckBox.isChecked():
            plotDialog.ui.plotCanvas.print_svg("element-Q-histograms.svg")
            self.ui.textBrowser.insertHtml("<H2>Per Element Atomic Charge Distributions<H2>\n")
            self.ui.textBrowser.insertHtml('<p><img src="element-Q-histograms.svg" /></p><br />\n')


    def handleBondLengthHistograms(self):
        plotDialog=BondLengthHistogramDialog(geometry=self.geo,parent=self)
        result=plotDialog.exec_()
        if result and plotDialog.ui.addOutputCheckBox.isChecked():
            plotDialog.ui.plotCanvas.print_svg("bond-length-histograms.svg")
            self.ui.textBrowser.insertHtml("<H2>Bond Length Distribution<H2>\n")
            self.ui.textBrowser.insertHtml('<p><img src="bond-length-histograms.svg" /></p><br />\n')
            
    
    def handleBondAngleHistograms(self):
        plotDialog=BondAngleHistogramDialog(geometry=self.geo,parent=self)
        result=plotDialog.exec_()
        if result and plotDialog.ui.addOutputCheckBox.isChecked():
            plotDialog.ui.plotCanvas.print_svg("bond-angle-histograms.svg")
            self.ui.textBrowser.insertHtml("<H2>Per Element Bond Angle Distributions<H2>\n")
            self.ui.textBrowser.insertHtml('<p><img src="bond-angle-histograms.svg" /></p><br />\n')


class BondAngleHistogramDialog(QtGui.QDialog):
    def __init__(self,geometry,parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = geostatspack4.plotDialog.Ui_Dialog()
        self.ui.setupUi(self)
        self.setWindowTitle("Bond Angle Histogram")
        
        self.geo=geometry
        
        self.ui.dataFileLineEdit.setText("bondangles.dat")

        self.ui.outputLayout=QtGui.QVBoxLayout(self.ui.outputGroupBox)
        self.ui.plotCanvas=MyMplCanvas(self.ui.outputGroupBox)
        self.ui.outputLayout.addWidget(self.ui.plotCanvas)
        
        self.ui.paramsLayout=QtGui.QHBoxLayout(self.ui.parametersGroupBox)
        self.ui.numberpars=QtGui.QFormLayout()
        self.ui.checkboxpars=QtGui.QVBoxLayout()
        self.ui.binCountLabel=QtGui.QLabel(self)
        self.ui.binCountLabel.setText("angle bins per center atom element")
        self.ui.binCountSpinBox=QtGui.QSpinBox(self)
        self.ui.binCountSpinBox.setValue(36)
        
        self.ui.numberpars.addRow(self.ui.binCountLabel, self.ui.binCountSpinBox)
        
        self.ui.spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.ui.checkboxpars.addItem(self.ui.spacerItem2)
        
        self.ui.paramsLayout.addLayout(self.ui.numberpars)
        self.ui.paramsLayout.addLayout(self.ui.checkboxpars)
        
        self.ui.outputGroupBox.layout=self.ui.outputLayout
        self.ui.parametersGroupBox.layout=self.ui.paramsLayout
        
        self.handlePlot()
        
        QtCore.QObject.connect(self.ui.saveButton, QtCore.SIGNAL("released()"),self.handleSave)
        QtCore.QObject.connect(self.ui.plotButton, QtCore.SIGNAL("released()"),self.handlePlot)
        
        
    def handlePlot(self):
        self.ui.plotCanvas.axes.cla()
        angleDict=self.geo.getBondAnglesByElementsHistograms(self.ui.binCountSpinBox.value())
        defaultColors=['r','g','b','c','m','k']
        colorcounter=0
        bonds=angleDict.keys()
        bonds.sort()
        for bond in bonds:
            label=self.geo.PTE[bond]
            width=angleDict[bond][1][1]-angleDict[bond][1][0]
            self.ui.plotCanvas.axes.bar(angleDict[bond][1][:-1],angleDict[bond][0],width,alpha=0.5,color=defaultColors[colorcounter%len(defaultColors)],label=label)
            colorcounter+=1
        self.ui.plotCanvas.axes.legend(loc=0)
        self.ui.plotCanvas.draw()
        
        

    def handleSave(self):
        self.handlePlot()
        if self.ui.saveDataCheckBox.isChecked():
            angleDict=self.geo.getBondAnglesByElementsHistograms(self.ui.binCountSpinBox.value())
            for bond in angleDict.keys():
                    label = self.geo.PTE[bond]
                    outfile=open(label+"."+str(self.ui.dataFileLineEdit.text()),"w")
                    print("#angle [deg]\tcount",file=outfile)
                    for i in range(len(angleDict[bond][0])): #@UnusedVariable
                        print("{0:e}\t{1:e}".format(angleDict[bond][1][i],angleDict[bond][0][i]),file=outfile)
                    outfile.close()
        self.accept()



class BondLengthHistogramDialog(QtGui.QDialog):
    def __init__(self,geometry,parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = geostatspack4.plotDialog.Ui_Dialog()
        self.ui.setupUi(self)
        self.setWindowTitle("Bond Length Histogram")
        
        self.geo=geometry
        
        self.ui.dataFileLineEdit.setText("bondlengths.dat")

        self.ui.outputLayout=QtGui.QVBoxLayout(self.ui.outputGroupBox)
        self.ui.plotCanvas=MyMplCanvas(self.ui.outputGroupBox)
        self.ui.outputLayout.addWidget(self.ui.plotCanvas)
        
        self.ui.paramsLayout=QtGui.QHBoxLayout(self.ui.parametersGroupBox)
        self.ui.numberpars=QtGui.QFormLayout()
        self.ui.checkboxpars=QtGui.QVBoxLayout()
        self.ui.binCountLabel=QtGui.QLabel(self)
        self.ui.binCountLabel.setText("length bins per bond type")
        self.ui.binCountSpinBox=QtGui.QSpinBox(self)
        self.ui.binCountSpinBox.setValue(30)
        self.ui.totalTooCheckBox=QtGui.QCheckBox(self)
        self.ui.totalTooCheckBox.setText("include total histogram")
        self.ui.totalTooCheckBox.setChecked(False)
        
        self.ui.numberpars.addRow(self.ui.binCountLabel, self.ui.binCountSpinBox)
        
        self.ui.checkboxpars.addWidget(self.ui.totalTooCheckBox)
        self.ui.spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.ui.checkboxpars.addItem(self.ui.spacerItem2)
        
        self.ui.paramsLayout.addLayout(self.ui.numberpars)
        self.ui.paramsLayout.addLayout(self.ui.checkboxpars)
        
        self.ui.outputGroupBox.layout=self.ui.outputLayout
        self.ui.parametersGroupBox.layout=self.ui.paramsLayout
        
        self.handlePlot()
        
        QtCore.QObject.connect(self.ui.saveButton, QtCore.SIGNAL("released()"),self.handleSave)
        QtCore.QObject.connect(self.ui.plotButton, QtCore.SIGNAL("released()"),self.handlePlot)
        
        
    def handlePlot(self):
        self.ui.plotCanvas.axes.cla()
        lengthDict=self.geo.getBondLengthHistograms(self.ui.binCountSpinBox.value())
        defaultColors=['r','g','b','c','m','k']
        colorcounter=0
        bonds=lengthDict.keys()
        bonds.sort()
        for bond in bonds:
            label=self.getBondLabel(bond)
            width=lengthDict[bond][0][1]-lengthDict[bond][0][0]
            if bond!=(-1,-1) or self.ui.totalTooCheckBox.isChecked()==True:
                self.ui.plotCanvas.axes.bar(lengthDict[bond][0],lengthDict[bond][1],width,alpha=0.5,color=defaultColors[colorcounter%len(defaultColors)],label=label)
                colorcounter+=1
        self.ui.plotCanvas.axes.legend(loc=0)
        self.ui.plotCanvas.draw()
        
        
        

    def getBondLabel(self, bond):
        if bond[0]==-1:
            return "total"
        else:
            return self.geo.PTE[bond[0]] + "-" + self.geo.PTE[bond[1]]


    def handleSave(self):
        self.handlePlot()
        if self.ui.saveDataCheckBox.isChecked():
            lengthDict=self.geo.getBondLengthHistograms(self.ui.binCountSpinBox.value())
            for bond in lengthDict.keys():
                if bond!=(-1,-1) or self.ui.totalTooCheckBox.isChecked()==True:
                    label = self.getBondLabel(bond)
                    outfile=open(label+"."+str(self.ui.dataFileLineEdit.text()),"w")
                    print("#r [a.u.]\tcount",file=outfile)
                    for i in range(len(lengthDict[bond][0])): #@UnusedVariable
                        print("{0:e}\t{1:e}".format(lengthDict[bond][0][i],lengthDict[bond][1][i]),file=outfile)
                    outfile.close()
        self.accept()



class ElementChargeHistogramDialog(QtGui.QDialog):
    def __init__(self,geometry,parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = geostatspack4.plotDialog.Ui_Dialog()
        self.ui.setupUi(self)
        self.setWindowTitle("Element Charge Histogram")
        
        self.geo=geometry
        
        self.ui.dataFileLineEdit.setText("elementcharges.dat")

        self.ui.outputLayout=QtGui.QVBoxLayout(self.ui.outputGroupBox)
        self.ui.plotCanvas=MyMplCanvas(self.ui.outputGroupBox)
        self.ui.outputLayout.addWidget(self.ui.plotCanvas)
        
        self.ui.paramsLayout=QtGui.QHBoxLayout(self.ui.parametersGroupBox)
        self.ui.numberpars=QtGui.QFormLayout()
        self.ui.checkboxpars=QtGui.QVBoxLayout()
        self.ui.binCountLabel=QtGui.QLabel(self)
        self.ui.binCountLabel.setText("charge bins per element")
        self.ui.binCountSpinBox=QtGui.QSpinBox(self)
        self.ui.binCountSpinBox.setValue(10)
        
        self.ui.numberpars.addRow(self.ui.binCountLabel, self.ui.binCountSpinBox)
        
        self.ui.spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.ui.checkboxpars.addItem(self.ui.spacerItem2)
        
        self.ui.paramsLayout.addLayout(self.ui.numberpars)
        self.ui.paramsLayout.addLayout(self.ui.checkboxpars)
        
        self.ui.outputGroupBox.layout=self.ui.outputLayout
        self.ui.parametersGroupBox.layout=self.ui.paramsLayout
        
        self.handlePlot()
        
        QtCore.QObject.connect(self.ui.saveButton, QtCore.SIGNAL("released()"),self.handleSave)
        QtCore.QObject.connect(self.ui.plotButton, QtCore.SIGNAL("released()"),self.handlePlot)
        
        
    def handlePlot(self):
        self.ui.plotCanvas.axes.cla()
        chargeDict=self.geo.elem_charges_hist(self.ui.binCountSpinBox.value())
        elementColors={
                       1:"#BEBEBE",
                       2:"b",
                       3:"c",
                       4:"m",
                       5:"brown",
                       6:"black",
                       7:"blue",
                       8:"r",
                       9:"g",
                       10:"b"}
        defaultColors=['r','g','b','c','m','k']
        for element in chargeDict.keys():
            width=chargeDict[element][0][1]-chargeDict[element][0][0]
            self.ui.plotCanvas.axes.bar(chargeDict[element][0],chargeDict[element][1],width,alpha=0.5,color=elementColors.get(element,defaultColors[element%len(defaultColors)]),label=self.geo.PTE[element])
        self.ui.plotCanvas.axes.legend(loc=0)
        self.ui.plotCanvas.draw()
        
        
        
    def handleSave(self):
        self.handlePlot()
        if self.ui.saveDataCheckBox.isChecked():
            chargeDict=self.geo.elem_charges_hist(self.ui.binCountSpinBox.value())
            for element in chargeDict.keys():
                    outfile=open(self.geo.PTE[element]+"."+str(self.ui.dataFileLineEdit.text()),"w")
                    print("#q [e]\tcount",file=outfile)
                    for i in range(len(chargeDict[element][0])): #@UnusedVariable
                        print("{0:e}\t{1:e}".format(chargeDict[element][0][i],chargeDict[element][1][i]),file=outfile)
                    outfile.close()
        self.accept()
        

    
class RDFDialog(QtGui.QDialog):
    def __init__(self,geometry,parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = geostatspack4.plotDialog.Ui_Dialog()
        self.ui.setupUi(self)
        self.setWindowTitle("Radial Distribution Function")
        
        self.geo=geometry
        
        self.ui.dataFileLineEdit.setText("rdf.dat")
        
        self.ui.outputLayout=QtGui.QVBoxLayout(self.ui.outputGroupBox)
        self.ui.plotCanvas=MyMplCanvas(self.ui.outputGroupBox)
        self.ui.outputLayout.addWidget(self.ui.plotCanvas)
        
        self.ui.paramsLayout=QtGui.QHBoxLayout(self.ui.parametersGroupBox)
        self.ui.numberpars=QtGui.QFormLayout()
        self.ui.checkboxpars=QtGui.QVBoxLayout()
        self.ui.binWidthLabel=QtGui.QLabel(self)
        self.ui.binWidthLabel.setText("step width (a.u.)")
        self.ui.bindWidthSpinBox=QtGui.QDoubleSpinBox(self)
        self.ui.bindWidthSpinBox.setValue(0.2)
        self.ui.bindWidthSpinBox.setSingleStep(0.05)
        self.ui.perElementCheckBox=QtGui.QCheckBox(self)
        self.ui.perElementCheckBox.setText("per element rdf")
        self.ui.perElementCheckBox.setChecked(False)
        self.ui.totalTooCheckBox=QtGui.QCheckBox(self)
        self.ui.totalTooCheckBox.setText("total rdf")
        self.ui.totalTooCheckBox.setChecked(True)
        self.ui.totalTooCheckBox.setEnabled(False)
        QtCore.QObject.connect(self.ui.perElementCheckBox, QtCore.SIGNAL("toggled(bool)"),
                               self.ui.totalTooCheckBox,QtCore.SLOT("setEnabled(bool)"))
        self.ui.spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        
        
        
        self.ui.numberpars.addRow(self.ui.binWidthLabel, self.ui.bindWidthSpinBox)
        
        self.ui.checkboxpars.addWidget(self.ui.perElementCheckBox)
        self.ui.checkboxpars.addWidget(self.ui.totalTooCheckBox)
        self.ui.checkboxpars.addItem(self.ui.spacerItem2)
        
        self.ui.paramsLayout.addLayout(self.ui.numberpars)
        self.ui.paramsLayout.addLayout(self.ui.checkboxpars)
        
        self.ui.outputGroupBox.layout=self.ui.outputLayout
        self.ui.parametersGroupBox.layout=self.ui.paramsLayout
        initRdf=self.geo.rdf()
        self.plotSingle(initRdf)
        
        QtCore.QObject.connect(self.ui.saveButton, QtCore.SIGNAL("released()"),self.handleSave)
        QtCore.QObject.connect(self.ui.plotButton, QtCore.SIGNAL("released()"),self.handlePlot)
        
 
        
    def plotSingle(self,rdf):
        """plot a single radial distribution function
        @type rdf: sequence of two sequences of scalars ([[r],[rdf]]) 
        @param rdf: sequences r and rdf(r)
        """
        self.ui.plotCanvas.axes.cla()
        self.ui.plotCanvas.axes.plot(rdf[0],rdf[1])
        self.ui.plotCanvas.draw()
    
    
    def plotMultiple(self,rdfs):
        """plot a single radial distribution function
        @type rdfs: sequence of dictionaries 
        @param rdfs: dictionaries containing at least "rdf" entries corresponding to plotSingle parameter, additionally can contain "label" string entry and "fill" boolean entry. 
        """
        self.ui.plotCanvas.axes.cla()
        plotargs=[]
        labels=[]
        for entry in rdfs:
            plotargs.append(entry["rdf"][0])
            plotargs.append(entry["rdf"][1])
            labels.append(entry.get("label","rdf"))
        self.ui.plotCanvas.axes.plot(*plotargs)
        self.ui.plotCanvas.axes.legend(labels)
        self.ui.plotCanvas.draw()
    
    
    def handleSave(self):
        self.handlePlot()
        if self.ui.saveDataCheckBox.isChecked():
            if self.ui.perElementCheckBox.isChecked():
                rdfs=self.getMultiRdfs()
                for entry in rdfs:
                    outfile=open(entry.get("label","rdf")+"."+str(self.ui.dataFileLineEdit.text()),"w")
                    print("#r [a.u.]\trdf(r)",file=outfile)
                    for i in range(len(entry["rdf"][0])):
                        print("{0:e}\t{1:e}".format(entry["rdf"][0][i],entry["rdf"][1][i]),file=outfile)
                    outfile.close()
            else:
                outfile=open(str(self.ui.dataFileLineEdit.text()),"w")
                print("#r [a.u.]\trdf(r)",file=outfile)
                rdf=self.geo.rdf(binwidth=self.ui.bindWidthSpinBox.value())
                for i in range(len(rdf[0])):
                    print("{0:e}\t{1:e}".format(rdf[0][i],rdf[1][i]),file=outfile)
                outfile.close()
        self.accept()
    
    

    def getMultiRdfs(self):
        rdfs = []
        if self.ui.totalTooCheckBox.isChecked():
            rdfs.append({"label":"total", "fillstyle":True, "rdf":self.geo.rdf(binwidth=self.ui.bindWidthSpinBox.value())})
        crossRdfs=self.geo.crossElementRDFs(self.ui.bindWidthSpinBox.value())
        for combo in crossRdfs.keys():
            lstring="{0:3>}-{1:3<}".format(self.geo.PTE[combo[0]],self.geo.PTE[combo[1]])
            rdfs.append({"label":lstring, "fillstyle":False, "rdf":crossRdfs[combo]})
        return rdfs

    def handlePlot(self):
        if self.ui.perElementCheckBox.isChecked():
            rdfs = self.getMultiRdfs()
            self.plotMultiple(rdfs)
        else:
            rdf=self.geo.rdf(binwidth=self.ui.bindWidthSpinBox.value())
            self.plotSingle(rdf)


class GeostatsHelp(QtGui.QDialog):
    def __init__(self,parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = geostatspack4.geostatsHelp.Ui_GeostatsHelp()
        self.ui.setupUi(self)


class PeriodicExpand(QtGui.QDialog):
    def __init__(self,parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = geostatspack4.periodicExpand.Ui_periodicExpand()
        self.ui.setupUi(self)


class chargeConstraintsDialog(QtGui.QDialog):
    def __init__(self,geo,parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = geostatspack4.dftbChargeConstraintDialog.Ui_Dialog()
        self.ui.setupUi(self)
        
        QtCore.QObject.connect(self.ui.constraintButton, QtCore.SIGNAL("released()"),
                               self.handleConstraint)
        
        self.geo=geo
        
    def handleConstraint(self):
        try:
            tokens=str(self.ui.serialListPlainTextEdit.toPlainText()).split()
            serials=[ int(iii)-1 for iii in tokens ]
        except ValueError as error:
            errorDialog=QtGui.QErrorMessage(parent=self)
            errorDialog.setModal(True)
            errorDialog.showMessage(error.args[0])
        except:
            errorDialog=QtGui.QErrorMessage(parent=self)
            errorDialog.setModal(True)
            errorDialog.showMessage("Error parsing list of atom serials")
            raise
        try:
            self.ui.constraintOutputTextBrowser.setPlainText(
                self.geo.getHSDChargeConstraintString(serials,self.ui.prefactorSpinBox.value()) )
        except comatsci.geometry.GeometryError as error:
            errorDialog=QtGui.QErrorMessage(parent=self)
            errorDialog.setModal(True)
            errorDialog.showMessage(error.args[0])
        except ValueError as error:
            errorDialog=QtGui.QErrorMessage(parent=self)
            errorDialog.setModal(True)
            errorDialog.showMessage(error.args[0])
            

class SimpleLinkDialog(QtGui.QDialog):
    def __init__(self, geo, parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = geostatspack4.embedDialogBase.Ui_Dialog()
        self.ui.setupUi(self)
        
        self.setWindowTitle("simple link atom generation")
        
        QtCore.QObject.connect(self.ui.buttonBox.button(QtGui.QDialogButtonBox.Apply), QtCore.SIGNAL("released()"),self.handleApply)
        QtCore.QObject.connect(self.ui.buttonBox, QtCore.SIGNAL("accepted()"), self.handleAccept)
        
        
        self.ui.neutralizeClusterCheckBox.setEnabled(False)
        
        self.geo=geo
        
        if len(self.geo.LayerDict.keys())<2:
            errorDialog=QtGui.QErrorMessage(parent=self)
            errorDialog.setModal(True)
            errorDialog.showMessage("Cannot link a Geometry with less than two layers!")
            self.reject()
        
        layerNames=QtCore.QStringList()
        for lIndex in geo.LayerDict.keys():
            layerNames.append(geo.LayerDict[lIndex].Name)
        
        self.ui.pchrComboBox.addItems(layerNames)
        pchrIndex=self.ui.pchrComboBox.findText("PCHR")
        if pchrIndex!=-1:
            self.ui.pchrComboBox.setCurrentIndex(pchrIndex)
        
        self.ui.qmzComboBox.addItems(layerNames)
        qmzIndex=self.ui.qmzComboBox.findText("QMZ")
        if qmzIndex!=-1:
            self.ui.qmzComboBox.setCurrentIndex(qmzIndex)
    
    
    def handleApply(self):
        QMZ=self.geo.layersubgeometry(self.geo.layerbyname(str(self.ui.qmzComboBox.currentText())))
        PCHR=self.geo.layersubgeometry(self.geo.layerbyname(str(self.ui.pchrComboBox.currentText())))
        self.linkedGeo,resultInfo=QMZ.SLALinkedGeometry(PCHR,distscale=self.ui.ladfSpinBox.value())
        self.ui.resultsTextBrowser.setPlainText(comatsci.utils.dictionaryPrettyPrint(resultInfo))
    
    
    def handleAccept(self):
        self.handleApply()
        self.accept()
        

class BCTCLinkDialog(QtGui.QDialog):
    def __init__(self, geo, parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = geostatspack4.embedDialogBase.Ui_Dialog()
        self.ui.setupUi(self)
        
        self.setWindowTitle("BCTC link atom generation")
        
        QtCore.QObject.connect(self.ui.buttonBox.button(QtGui.QDialogButtonBox.Apply), QtCore.SIGNAL("released()"),self.handleApply)
        QtCore.QObject.connect(self.ui.buttonBox, QtCore.SIGNAL("accepted()"), self.handleAccept)
        
        
        self.ui.neutralizeClusterCheckBox.setEnabled(True)
        
        self.ui.bctcCoeffLayoutWidget=QtGui.QWidget(self)
        self.ui.bctcCoeffLayout = QtGui.QHBoxLayout(self.ui.bctcCoeffLayoutWidget)
        self.ui.bctcCoeffLayout.setObjectName("bctcCoeffLayout")
        
        self.ui.bctcCoeffLabel=QtGui.QLabel("BCTC cond charge transfer coefficients",self.ui.verticalLayoutWidget)
        self.ui.verticalLayout.insertWidget(1,self.ui.bctcCoeffLabel)
        
        
        
        self.ui.bctcCoeffButtonLayoutWidget = QtGui.QWidget(self)
        self.ui.bctcCoeffButtonLayout = QtGui.QVBoxLayout(self.ui.bctcCoeffButtonLayoutWidget)
        self.ui.bctcCoeffButtonLayout.setObjectName("bctcCoeffButtonLayout")
        
        self.ui.bctcCoeffFileButton=QtGui.QPushButton("read file",self.ui.bctcCoeffButtonLayoutWidget)
        self.ui.bctcCoeffCalcButton=QtGui.QPushButton("calculate",self.ui.bctcCoeffButtonLayoutWidget)
        self.ui.bctcCoeffButtonLayout.addWidget(self.ui.bctcCoeffFileButton)
        self.ui.bctcCoeffButtonLayout.addWidget(self.ui.bctcCoeffCalcButton)
        
        QtCore.QObject.connect(self.ui.bctcCoeffCalcButton, QtCore.SIGNAL("released()"), self.handleCalculate)
        QtCore.QObject.connect(self.ui.bctcCoeffFileButton, QtCore.SIGNAL("released()"), self.handleReadFile)
        
        self.ui.bctcCoeffTable=QtGui.QTableWidget(self)
        
        self.ui.bctcCoeffLayout.insertWidget(-1,self.ui.bctcCoeffButtonLayoutWidget)
        self.ui.bctcCoeffLayout.addWidget(self.ui.bctcCoeffTable)
        
        
        self.ui.verticalLayout.insertWidget(2,self.ui.bctcCoeffLayoutWidget)
        self.ui.verticalLayout.removeItem(self.ui.verticalLayout.findChild(QtGui.QSpacerItem, QtCore.QString("spacerItem")))
        
        self.geo=geo
        self.chargeTransfers=None
        
        symlist,symdict=self.geo.getatomsymlistdict()
        symlist.sort()
        
        self.ui.bctcCoeffTable.setColumnCount(2)
        self.ui.bctcCoeffTable.setHorizontalHeaderLabels(QtCore.QStringList(["elements","coefficient"]))
        
        del symdict
        # now construct list of _unique_ combinations and a dictionary mapping element-element pairs to combination list indices
        elementCombinations=[]
        self.NameElementCombinations={}
        for i in range(len(symlist)):
            for j in range (i,len(symlist)):
                elementCombinations.append((symlist[i],symlist[j]))
                self.NameElementCombinations["{0:s}-{1:s}".format(self.geo.PTE[symlist[i]],self.geo.PTE[symlist[j]])]=(symlist[i],symlist[j])
        
        self.ui.bctcCoeffTable.setRowCount(len(self.NameElementCombinations.keys()))
        
        for i in range(len(self.NameElementCombinations.keys())):
            nameItem=QtGui.QTableWidgetItem()
            nameItem.setText(self.NameElementCombinations.keys()[i])
            valueItem=QtGui.QTableWidgetItem()
            valueItem.setText("0.0")
            self.ui.bctcCoeffTable.setItem(i,0,nameItem)
            self.ui.bctcCoeffTable.setItem(i,1,valueItem)
            nameItem.setFlags(nameItem.flags() & ~QtCore.Qt.ItemIsEditable)
            
        
        QtCore.QObject.connect(self.ui.bctcCoeffTable, QtCore.SIGNAL("cellChanged(int,int)"), self.tableToChargeTransfers)
        
        
        if len(self.geo.LayerDict.keys())<2:
            errorDialog=QtGui.QErrorMessage(parent=self)
            errorDialog.setModal(True)
            errorDialog.showMessage("Cannot link a Geometry with less than two layers!")
            self.reject()
        
        layerNames=QtCore.QStringList()
        for lIndex in geo.LayerDict.keys():
            layerNames.append(geo.LayerDict[lIndex].Name)
        
        self.ui.pchrComboBox.addItems(layerNames)
        pchrIndex=self.ui.pchrComboBox.findText("PCHR")
        if pchrIndex!=-1:
            self.ui.pchrComboBox.setCurrentIndex(pchrIndex)
        
        self.ui.qmzComboBox.addItems(layerNames)
        qmzIndex=self.ui.qmzComboBox.findText("QMZ")
        if qmzIndex!=-1:
            self.ui.qmzComboBox.setCurrentIndex(qmzIndex)
    
    
    def handleApply(self):
        QMZ=self.geo.layersubgeometry(self.geo.layerbyname(str(self.ui.qmzComboBox.currentText())))
        PCHR=self.geo.layersubgeometry(self.geo.layerbyname(str(self.ui.pchrComboBox.currentText())))
        self.linkedGeo,resultInfo=QMZ.BCTCLinkedGeometry(PCHR,distscale=self.ui.ladfSpinBox.value(),chargeTransfers=self.chargeTransfers,neutralize=self.ui.neutralizeClusterCheckBox.isChecked())
        self.ui.resultsTextBrowser.setPlainText(comatsci.utils.dictionaryPrettyPrint(resultInfo))
    
    

    def chargeTransfersToTable(self, tempCTC):
        for row in range(self.ui.bctcCoeffTable.rowCount()):
            rowName = str(self.ui.bctcCoeffTable.item(row, 0).text())
            self.ui.bctcCoeffTable.item(row, 1).setText("{0:10.6f}".format(tempCTC[self.NameElementCombinations[rowName]]))
            

    def handleReadFile(self):
        bctcFileName=str(QtGui.QFileDialog.getOpenFileName(parent=self, caption=QtCore.QString("Open BCTC Coefficients"),
                                                      filter=QtCore.QString("BCTC coefficient files (*.bctc);;all supported files (*.*)")))
        if bctcFileName!="":
            tempCoefficients={}
            bctcFile=open(str(bctcFileName),"r")
            for line in list(bctcFile):
                if not line.strip()[0] in (";","#"):
                    tokens=line.split()
                    try:
                        elementA=int(tokens[0])
                        elementB=int(tokens[1])
                        coefficient=float(tokens[2])
                    except ValueError as error:
                        errorDialog=QtGui.QErrorMessage(parent=self)
                        errorDialog.setModal(True)
                        errorDialog.showMessage("Error parsing file {0:s}: {1:s}".format(str(bctcFileName),error.args[0]))
                        raise
                    except:
                        errorDialog=QtGui.QErrorMessage(parent=self)
                        errorDialog.setModal(True)
                        errorDialog.showMessage("Error parsing file {0:s}".format(str(bctcFileName)))
                        raise
                    tempCoefficients[(elementA,elementB)]=coefficient
            bctcFile.close()
            self.chargeTransfersToTable(tempCoefficients)
            self.chargeTransfers=tempCoefficients


    def handleCalculate(self):
        tempCTC=self.geo.getElementElementChargeTransfers()[0]
        self.chargeTransfersToTable(tempCTC)
        self.chargeTransfers=tempCTC
        
    
    
    def handleAccept(self):
        self.tableToChargeTransfers()
        self.handleApply()
        self.accept()
        
    
    def tableToChargeTransfers(self):
        oldChargeTransfers=self.chargeTransfers
        self.chargeTransfers={}
        for row in range(self.ui.bctcCoeffTable.rowCount()):
            rowName=str(self.ui.bctcCoeffTable.item(row, 0).text())
            try:
                rowValue=float(str(self.ui.bctcCoeffTable.item(row, 1).text()))
            except ValueError as error:
                errorDialog=QtGui.QErrorMessage(parent=self)
                errorDialog.setModal(True)
                errorDialog.showMessage(error.args[0])
                self.chargeTransfersToTable(oldChargeTransfers)
                raise
            except:
                errorDialog=QtGui.QErrorMessage(parent=self)
                errorDialog.setModal(True)
                errorDialog.showMessage("Error parsing Charge Coefficient Row")
                raise
            self.chargeTransfers[self.NameElementCombinations[rowName]]=rowValue
        
        
class EditAtomsDialog(QtGui.QDialog):
    def __init__(self, geo, parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = geostatspack4.editAtoms.Ui_Dialog()
        self.ui.setupUi(self)
        
        self.geo=copy.deepcopy(geo)
        self.originalGeo=geo
        
        self.gTM=GeometryTableModel(self.geo)
        
        
        self.ui.atomsTableView.setModel(self.gTM)
        self.ui.atomsTableView.setItemDelegate(GeometryAtomDelegate(self))
        
        QtCore.QObject.connect(self.ui.addAtomButton, QtCore.SIGNAL("released()"),self.handleAddAtom)
        QtCore.QObject.connect(self.ui.delAtomButton, QtCore.SIGNAL("released()"),self.handleDeleteAtom)
        QtCore.QObject.connect(self.ui.buttonBox.button(QtGui.QDialogButtonBox.Reset), QtCore.SIGNAL("released()"),self.handleReset)
        
        self.ui.atomsTableView.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        
    
    
    def handleReset(self):
        self.geo=copy.deepcopy(self.originalGeo)
        self.gTM.geo=self.geo
        self.gTM.reset()
        
    
    def handleAddAtom(self):
        self.geo.addatom(0,(0.0,0.0,0.0))
        self.gTM.reset()
        self.ui.atomsTableView.scrollToBottom()
    
    
    def handleDeleteAtom(self):
        atomIndex=self.ui.atomsTableView.currentIndex().row()
        if atomIndex >= 0 and atomIndex < self.geo.Atomcount:
            self.geo.delatom(atomIndex)
        self.gTM.reset()
    
        
    
class GeometryTableModel(QtCore.QAbstractTableModel):
    """
    Qt Model of Molecular Geometries
    """
    
    ELEMENT=0
    SUBTYPE=1
    X=2
    Y=3
    Z=4
    LAYER=5
    CHARGE=6
    
    
    def __init__(self,geo):
        """Construct Model for Geometry
        @type geo: comatsci.Geometry
        @param geo: the molecular geometry object to represent
        """
        super(GeometryTableModel,self).__init__()
        self.geo=geo
        self.dirty=False
        
    
    def rowCount(self, index=QtCore.QModelIndex()):
        return self.geo.Atomcount
    
    
    def columnCount(self, index=QtCore.QModelIndex()):
        return 7
    
    
    def data(self, index, role=QtCore.Qt.DisplayRole):
        if not index.isValid() or not (0 <= index.row() < self.geo.Atomcount):
            return QtCore.QVariant()
        row=index.row()
        column=index.column()
        if role==QtCore.Qt.DisplayRole:
            if column==self.ELEMENT:
                return QtCore.QVariant(comatsci.geometry.Geometry.PTE[self.geo.AtomTypes[row]])
            elif column==self.SUBTYPE:
                return QtCore.QVariant(self.geo.AtomSubTypes[row])
            elif column==self.X:
                return QtCore.QVariant(float(self.geo.Geometry[row][0]*comatsci.constants.ANGSTROM))
            elif column==self.Y:
                return QtCore.QVariant(float(self.geo.Geometry[row][1]*comatsci.constants.ANGSTROM))
            elif column==self.Z:
                return QtCore.QVariant(float(self.geo.Geometry[row][2]*comatsci.constants.ANGSTROM))
            elif column==self.LAYER:
                return QtCore.QVariant(self.geo.LayerDict[self.geo.AtomLayers[row]].Name)
            elif column==self.CHARGE:
                return QtCore.QVariant(self.geo.AtomCharges[row])
        elif role==QtCore.Qt.TextAlignmentRole:
            return QtCore.QVariant(int(QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter))
        elif role==QtCore.Qt.TextColorRole:
            if column==self.CHARGE:
                if self.geo.AtomCharges[row]==0:
                    return QtCore.QVariant(QtGui.QColor(QtCore.Qt.black))
                if self.geo.AtomCharges[row]<0:
                    return QtCore.QVariant(QtGui.QColor(QtCore.Qt.red))
                if self.geo.AtomCharges[row]>0:
                    return QtCore.QVariant(QtGui.QColor(QtCore.Qt.blue))
        # chatch-all default return
        return QtCore.QVariant()
    
    
    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role==QtCore.Qt.TextAlignmentRole:
            if orientation==QtCore.Qt.Horizontal:
                return QtCore.QVariant(int(QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter))
            else:
                return QtCore.QVariant(int(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter))
        elif role!= QtCore.Qt.DisplayRole:
            return QtCore.QVariant()
        if orientation==QtCore.Qt.Horizontal:
            if section==self.ELEMENT:
                return QtCore.QVariant("Element")
            elif section==self.SUBTYPE:
                return QtCore.QVariant("Type")
            elif section==self.X:
                return QtCore.QVariant("X [A]")
            elif section==self.Y:
                return QtCore.QVariant("Y [A]")
            elif section==self.Z:
                return QtCore.QVariant("Z [A]")
            elif section==self.LAYER:
                return QtCore.QVariant("Layer")
            elif section==self.CHARGE:
                return QtCore.QVariant("Charge [e]")
            else:
                return QtCore.QVariant(int(section+1))
                
                
                
    def flags(self,index):
        if not index.isValid():
            return QtCore.Qt.ItemIsEnabled
        elif index.column() in (self.ELEMENT,self.SUBTYPE,self.X,self.Y,self.Z,self.CHARGE,self.LAYER):
            return QtCore.Qt.ItemFlags(QtCore.QAbstractTableModel.flags(self,index)|QtCore.Qt.ItemIsEditable)
        else:
            return QtCore.Qt.ItemFlags(QtCore.QAbstractTableModel.flags(self,index))
        
        
        
        
    def setData(self, index, value, role=QtCore.Qt.EditRole):
        if index.isValid() and 0 < index.row() < self.geo.Atomcount:
            atom=index.row()
            column=index.column()
            if column==self.ELEMENT:
                self.geo.AtomTypes[atom]=self.geo.RPTE[str(value.toString()).lower()]
            if column==self.SUBTYPE:
                self.geo.AtomSubTypes[atom]=value.toString()
            if column==self.X:
                self.geo.Geometry[atom][0]=float(value.toString())/comatsci.constants.ANGSTROM
            if column==self.Y:
                self.geo.Geometry[atom][1]=float(value.toString())/comatsci.constants.ANGSTROM
            if column==self.Z:
                self.geo.Geometry[atom][2]=float(value.toString())/comatsci.constants.ANGSTROM
            if column==self.CHARGE:
                self.geo.AtomCharges[atom]=float(value.toString())
            if column==self.LAYER:
                layerName=str(value.toString())
                if layerName in self.geo.layerNames():
                    self.geo.AtomLayers[atom]=self.geo.layerbyname(layerName)
                else:
                    newlayer=self.geo.addlayer(layerName)
                    self.geo.AtomLayers[atom]=newlayer
            self.dirty=True
            self.emit(QtCore.SIGNAL("dataChanged(QModelIndex.QModelIndex)"),index,index)
            return True
        return False
    
    
class GeometryAtomDelegate(QtGui.QItemDelegate):
    
    def createEditor(self, parent, option, index):
        if index.column()==GeometryTableModel.ELEMENT:
            combobox=QtGui.QComboBox(parent)
            combobox.addItems(comatsci.geometry.Geometry.PTE)
            return combobox
        elif index.column()==GeometryTableModel.SUBTYPE:
            combobox=QtGui.QComboBox(parent)
            combobox.addItems(QtCore.QStringList(dict.fromkeys(index.model().geo.AtomSubTypes).keys()))
            combobox.setEditable(True)
            return combobox
        elif index.column()==GeometryTableModel.LAYER:
            combobox=QtGui.QComboBox(parent)
            combobox.addItems(index.model().geo.layerNames())
            print(index.model().geo.layerNames())
            combobox.setEditable(True)
            return combobox
        if index.column() in (GeometryTableModel.X,GeometryTableModel.Y,GeometryTableModel.Z):
            spinbox=QtGui.QDoubleSpinBox(parent)
            spinbox.setSingleStep(0.1)
            spinbox.setMinimum(-100000.0)
            spinbox.setMaximum(+100000.0)
            spinbox.setDecimals(6)
            return spinbox
        if index.column()==GeometryTableModel.CHARGE:
            spinbox=QtGui.QDoubleSpinBox(parent)
            spinbox.setSingleStep(0.01)
            spinbox.setMinimum(-100.0)
            spinbox.setMaximum(+100.0)
            spinbox.setDecimals(6)
            return spinbox
        else:
            return QtGui.QItemDelegate.createEditor(self, parent, option, index)
        
        
    def setEditorData(self, editor, index):
        text=index.model().data(index, QtCore.Qt.DisplayRole).toString()
        if index.column()==GeometryTableModel.ELEMENT:
            editor.setCurrentIndex(editor.findText(text))
        elif index.column==GeometryTableModel.SUBTYPE:
            editor.setCurrentIndex(editor.findText(text))
        if index.column() in (GeometryTableModel.X,GeometryTableModel.Y,GeometryTableModel.Z,GeometryTableModel.CHARGE):
            editor.setValue(float(text))
        else:
            return QtGui.QItemDelegate.setEditorData(self, editor, index)
        
        
    
    def setModelData(self, editor, model, index):
        if index.column() in (GeometryTableModel.X,GeometryTableModel.Y,GeometryTableModel.Z,GeometryTableModel.CHARGE):
            model.setData(index,QtCore.QVariant(editor.value()))
        elif index.column() in (GeometryTableModel.ELEMENT, GeometryTableModel.SUBTYPE, GeometryTableModel.LAYER):
            model.setData(index,QtCore.QVariant(editor.currentText()))


def mainfunc():
    global app, argList, geostatsWidget
    app = QtGui.QApplication(sys.argv)
    argList = app.arguments()
    geostatsWidget = MainWindow()
    print(list(argList))
    # if argList.count()==2:
    try:
        geostatsWidget.readGeometry(str(argList.last()))
    except comatsci.geometry.GeometryError:
        pass
    #elif argList.count()>2:
    #    raise ValueError("More than one geometry file specified on command line. Abort.")
    geostatsWidget.show()
    sys.exit(app.exec_())


if __name__=="__main__":
    mainfunc()
    