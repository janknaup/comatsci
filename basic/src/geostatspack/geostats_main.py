# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'geostats_main.ui'
#
# Created: Do Jul 16 16:28:05 2009
#      by: The PyQt User Interface Compiler (pyuic) 3.17.6
#
# WARNING! All changes made in this file will be lost!


from qt import *
##############################################################################
# geostats_main.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup , Knaup@bccms.uni-bremen.de
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

import comatsci
import EditGeometry
import periodicExpand
import elementChargeHistograms
import rdfWizard
import geostatshelp
import linkLayersSLA
import linkLayersHCS
import linkLayersBCTC
import chargeConstraintsWizard

image0_data = \
    "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d" \
    "\x49\x48\x44\x52\x00\x00\x00\x14\x00\x00\x00\x14" \
    "\x08\x06\x00\x00\x00\x8d\x89\x1d\x0d\x00\x00\x02" \
    "\x0e\x49\x44\x41\x54\x38\x8d\xad\x94\xbf\x4b\x5b" \
    "\x51\x14\xc7\x3f\x37\xbe\x67\x22\xf1\xbd\xa4\xa5" \
    "\x04\x89\x10\x77\x71\x2e\xe4\x2f\x70\x34\xb8\x5a" \
    "\x87\x2e\x2e\xfe\x01\x0e\x6e\xce\x4a\xa6\xd0\x35" \
    "\xae\x01\xc1\xb5\xb4\x5d\xb2\x26\x10\x44\x0d\x91" \
    "\x2c\x31\xf5\xa5\xf2\xe8\x7b\x49\xa3\xc1\x34\x2f" \
    "\xf6\x74\x89\xf9\x61\x94\xda\xd7\x7e\xe1\xc2\xbd" \
    "\xe7\xc2\xe7\x9c\xfb\x3d\x87\xab\x44\x04\x80\x7a" \
    "\xbd\x2e\x37\x37\x37\x18\x86\x41\x2c\x16\x23\x14" \
    "\x0a\x29\x7c\x28\x00\xa0\x94\x92\xe3\xe3\x63\xba" \
    "\xdd\x2e\xba\xae\x53\x2c\x16\x91\x87\x4c\x7e\x80" \
    "\x00\xb6\x6d\x63\xdb\x36\xa5\x52\x89\xd3\xd3\x53" \
    "\x9a\xcd\xa6\x1f\x1e\xda\xc3\xc6\xb2\x2c\xf2\xf9" \
    "\x3c\xbd\x5e\x0f\xd7\x75\xd9\xda\xda\xf2\x05\x0c" \
    "\x28\xa5\x04\x20\x95\x4a\xe1\x38\x0e\xb6\x6d\xb3" \
    "\xb1\xb1\xc1\xcc\xcc\x8c\x2f\x0f\x95\x88\xa0\x94" \
    "\x1a\x5a\x56\x2e\x97\xa9\x56\xab\x24\x93\x49\x62" \
    "\xb1\xd8\x5f\x43\x03\x8f\x03\xcb\xcb\xcb\x74\x2e" \
    "\x2f\xb1\x36\x37\xe9\xd4\xeb\xd2\x3e\x3f\x97\xfb" \
    "\x6e\xf7\xe5\x0d\x1a\x54\x26\xe3\xba\xaa\x54\xe4" \
    "\x83\xae\xcb\x67\x90\x33\xa5\xe4\x2c\x1a\x95\x7a" \
    "\x26\x23\x22\xc2\x9f\x96\xf6\x54\x92\x8f\xfb\xfb" \
    "\x04\x3d\x8f\x30\x30\xaf\x69\xbc\x0a\x87\x29\x6f" \
    "\x6f\xf3\x6b\x71\x51\x96\xd6\xd6\x86\x36\x34\x1a" \
    "\x0d\xb9\xb8\xb8\x60\x61\x61\x81\x50\x28\x44\x22" \
    "\x91\x18\x01\x4d\xd7\x04\xe0\x3a\xf8\x8d\x2f\x87" \
    "\x87\xbc\x05\xda\x9a\xc6\xfc\xec\x2c\xdf\x2d\x8b" \
    "\x0e\x50\xd8\xdb\xc3\xbc\xbb\x13\x56\x57\x01\xc8" \
    "\x66\xb3\xac\xac\xac\x10\x8f\xc7\x69\xb5\x5a\xdc" \
    "\xde\xde\x4e\x57\xd8\x6e\x36\xb9\xf7\x3c\x1c\xe0" \
    "\xab\xa6\xc1\xdc\x1c\xcd\x4e\x87\x36\x50\x3d\x39" \
    "\xc1\xda\xd9\x41\xae\xaf\x01\xa8\xd5\x6a\x04\x83" \
    "\x41\xfa\xfd\x3e\xae\xeb\x12\x89\x44\x50\x80\x00" \
    "\x18\x8e\x01\xc0\x8f\x68\x8b\x77\x4b\x4b\xbc\xb9" \
    "\xba\x22\x32\xe8\xda\x4f\x20\x33\xb8\x6f\xbf\x6e" \
    "\x0f\x9f\x9c\x4e\xa7\xa5\x52\xa9\xa0\xeb\x3a\x96" \
    "\x65\xb1\xbb\xbb\x3b\x6a\x8a\xe1\x18\x62\x38\x86" \
    "\x88\x88\x7c\x3a\x3a\x1a\x9e\xdf\x83\xac\x4f\xde" \
    "\x0f\x1b\xe0\x79\x9e\xe4\x72\x39\x39\x38\x38\x90" \
    "\x42\xa1\x20\x22\x32\x9a\x43\x63\x54\xc1\x84\xa7" \
    "\x53\x96\x8c\x55\xf8\x94\xa6\x80\x2f\xd5\x73\xe0" \
    "\xa9\xc1\x7e\xa9\x4c\xd7\x14\xd3\x35\xa7\x06\xde" \
    "\x37\xf0\x39\xf0\x3f\x03\x1f\x6b\x38\x36\xff\x43" \
    "\x22\xa2\xc6\x7f\x9b\x09\x93\xfd\xc6\x7e\x03\x76" \
    "\x35\x47\x87\x2d\xbf\x85\xbe\x00\x00\x00\x00\x49" \
    "\x45\x4e\x44\xae\x42\x60\x82"
image1_data = \
    "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d" \
    "\x49\x48\x44\x52\x00\x00\x00\x16\x00\x00\x00\x16" \
    "\x08\x06\x00\x00\x00\xc4\xb4\x6c\x3b\x00\x00\x00" \
    "\x99\x49\x44\x41\x54\x38\x8d\xed\x94\x41\x0e\x85" \
    "\x20\x0c\x44\x5f\x89\xc7\x36\x7f\x61\xbc\x77\x5d" \
    "\x28\x48\xa4\x28\x60\xff\xce\xd9\x54\x8b\xbe\x8e" \
    "\x13\x04\x3e\x1d\x92\x81\x77\xf4\x81\xa1\x23\xdc" \
    "\x2b\x34\xf6\xf4\x7a\x3d\xe2\xb8\x65\xa8\x84\x3f" \
    "\x40\x01\x98\x2a\x0b\x3d\x5f\x62\xc5\x83\x00\xaa" \
    "\x1a\xd7\x05\x50\x44\x9a\xb9\xd5\x07\xa7\x73\xa8" \
    "\xa4\xba\x4f\x92\xa2\xdf\x33\x3c\x64\xc6\x3b\xeb" \
    "\xbd\x82\xe5\xb8\xad\xde\xcb\xcc\x78\x20\xeb\x42" \
    "\x66\xc6\x39\x74\x5d\xfa\x80\xf3\x6f\xaf\x66\xc6" \
    "\x6f\xa1\x9c\x3f\x88\x2f\xb4\x70\xec\x05\xcd\xc0" \
    "\xbe\xd0\x78\x93\xf6\x8e\x17\x14\x92\x63\x5f\x68" \
    "\x6c\x3e\xef\xf6\xba\x3c\x8f\xdd\x36\x6d\xc4\xc0" \
    "\x45\x2c\x87\x81\xf8\x08\x00\x00\x00\x00\x49\x45" \
    "\x4e\x44\xae\x42\x60\x82"
image2_data = \
    "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d" \
    "\x49\x48\x44\x52\x00\x00\x00\x16\x00\x00\x00\x16" \
    "\x08\x06\x00\x00\x00\xc4\xb4\x6c\x3b\x00\x00\x00" \
    "\xa0\x49\x44\x41\x54\x38\x8d\xd5\x95\x4d\x0a\x80" \
    "\x20\x10\x85\x9f\xd1\x46\x68\xe1\x8d\xe6\x62\xd2" \
    "\x22\xbc\x98\x37\x6a\x21\xb4\xac\x45\x19\x92\xc6" \
    "\x64\x69\xe0\xb7\xf1\x87\xf1\xf1\x1c\x47\x05\x2a" \
    "\x21\x8e\x76\x2d\xad\xdb\xfb\x9e\x99\xf6\x56\x8f" \
    "\x80\xb5\x36\x4b\x85\x88\xce\x35\x44\x04\x00\xe8" \
    "\x0a\x39\x8c\xe8\xf9\x90\x34\xd2\x29\x2c\xc3\x7c" \
    "\x8e\xbd\x53\x0f\xeb\x58\x3a\x05\xe9\x54\x34\x1f" \
    "\x8a\x02\x7b\x2a\x7d\x3a\x1f\x09\xbf\x85\x4d\xc5" \
    "\xd5\xd9\x53\xaa\x39\x6e\x4f\x38\xca\xb1\x99\xe2" \
    "\xd2\xe1\x08\xab\xe1\x56\xf8\x2e\x30\x97\x7f\xcb" \
    "\x4d\x8f\xf9\x42\xd7\x5d\xbe\xbe\xd2\xe1\x43\x95" \
    "\x3a\x93\xf6\xca\xad\x3d\x61\x11\xf4\x4b\x7d\x4f" \
    "\x82\x0f\xf9\xc0\x06\x9b\xb5\x1e\xcd\xed\x31\x8c" \
    "\x5c\x00\x00\x00\x00\x49\x45\x4e\x44\xae\x42\x60" \
    "\x82"

class MainWindow(QMainWindow):
    def __init__(self,parent = None,name = None,fl = 0):
        QMainWindow.__init__(self,parent,name,fl)
        self.statusBar()

        self.image0 = QPixmap()
        self.image0.loadFromData(image0_data,"PNG")
        self.image1 = QPixmap()
        self.image1.loadFromData(image1_data,"PNG")
        self.image2 = QPixmap()
        self.image2.loadFromData(image2_data,"PNG")
        if not name:
            self.setName("geostats")

        self.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed,0,0,self.sizePolicy().hasHeightForWidth()))
        self.setMinimumSize(QSize(21,169))
        self.setMaximumSize(QSize(530,440))
        self.setBaseSize(QSize(530,440))
        self.setIcon(self.image0)

        self.setCentralWidget(QWidget(self,"qt_central_widget"))

        self.frame3 = QFrame(self.centralWidget(),"frame3")
        self.frame3.setGeometry(QRect(0,340,530,50))
        self.frame3.setFrameShape(QFrame.StyledPanel)
        self.frame3.setFrameShadow(QFrame.Raised)

        self.periodicExpandButton = QPushButton(self.frame3,"periodicExpandButton")
        self.periodicExpandButton.setEnabled(0)
        self.periodicExpandButton.setGeometry(QRect(140,10,120,31))

        self.editAtomsButton = QPushButton(self.frame3,"editAtomsButton")
        self.editAtomsButton.setEnabled(0)
        self.editAtomsButton.setGeometry(QRect(10,10,120,31))

        self.VMDButton = QPushButton(self.frame3,"VMDButton")
        self.VMDButton.setEnabled(0)
        self.VMDButton.setGeometry(QRect(399,10,120,31))

        self.whateverButton = QPushButton(self.frame3,"whateverButton")
        self.whateverButton.setEnabled(0)
        self.whateverButton.setGeometry(QRect(270,10,120,31))

        self.textBrowser1 = QTextBrowser(self.centralWidget(),"textBrowser1")
        self.textBrowser1.setGeometry(QRect(0,0,530,330))
        textBrowser1_font = QFont(self.textBrowser1.font())
        textBrowser1_font.setFamily("Monospace")
        self.textBrowser1.setFont(textBrowser1_font)

        self.fileOpenAction = QAction(self,"fileOpenAction")
        self.fileOpenAction.setIconSet(QIconSet(self.image1))
        self.fileSaveAction = QAction(self,"fileSaveAction")
        self.fileSaveAction.setEnabled(0)
        self.fileSaveAction.setIconSet(QIconSet(self.image2))
        self.fileSaveAsAction = QAction(self,"fileSaveAsAction")
        self.fileSaveAsAction.setEnabled(0)
        self.fileExitAction = QAction(self,"fileExitAction")
        self.helpContentsAction = QAction(self,"helpContentsAction")
        self.helpAboutAction = QAction(self,"helpAboutAction")
        self.statisticsCoordinationsAction = QAction(self,"statisticsCoordinationsAction")
        self.statisticsCoordinationsAction.setEnabled(0)
        self.statisticsSave_StatisticsAction = QAction(self,"statisticsSave_StatisticsAction")
        self.statisticsSave_StatisticsAction.setEnabled(0)
        self.statisticsbond_listAction = QAction(self,"statisticsbond_listAction")
        self.statisticsbond_listAction.setEnabled(0)
        self.statisticsallAction = QAction(self,"statisticsallAction")
        self.editedit_atomsAction = QAction(self,"editedit_atomsAction")
        self.editedit_atomsAction.setEnabled(0)
        self.editperiodic_expandAction = QAction(self,"editperiodic_expandAction")
        self.editperiodic_expandAction.setEnabled(0)
        self.fileView_in_VMDAction = QAction(self,"fileView_in_VMDAction")
        self.fileView_in_VMDAction.setEnabled(0)
        self.statisticsChargesAction = QAction(self,"statisticsChargesAction")
        self.statisticsChargesAction.setEnabled(0)
        self.editgenerate_linkatomsAction = QAction(self,"editgenerate_linkatomsAction")
        self.editgenerate_linkatomsAction.setEnabled(0)
        self.statisticsgraphsAction = QAction(self,"statisticsgraphsAction")
        self.statisticsgraphsAction.setEnabled(0)
        self.statisticsgraphselement_charge_histogramsAction = QAction(self,"statisticsgraphselement_charge_histogramsAction")
        self.statisticsgraphselement_charge_histogramsAction.setEnabled(0)
        self.statisticsgraphsradial_distribution_functionAction = QAction(self,"statisticsgraphsradial_distribution_functionAction")
        self.statisticsgraphsradial_distribution_functionAction.setEnabled(0)
        self.statisticssave_BCT_coefficientsAction = QAction(self,"statisticssave_BCT_coefficientsAction")
        self.statisticssave_BCT_coefficientsAction.setEnabled(0)
        self.statisticsget_charge_constraintsAction = QAction(self,"statisticsget_charge_constraintsAction")
        self.statisticsget_charge_constraintsAction.setEnabled(0)
        self.editBCTC_link_LayersAction = QAction(self,"editBCTC_link_LayersAction")
        self.editBCTC_link_LayersAction.setEnabled(0)
        self.editHCS_link_layersAction = QAction(self,"editHCS_link_layersAction")
        self.editHCS_link_layersAction.setEnabled(0)




        self.MenuBar = QMenuBar(self,"MenuBar")


        self.fileMenu = QPopupMenu(self)
        self.fileOpenAction.addTo(self.fileMenu)
        self.fileSaveAction.addTo(self.fileMenu)
        self.fileSaveAsAction.addTo(self.fileMenu)
        self.fileMenu.insertSeparator()
        self.fileView_in_VMDAction.addTo(self.fileMenu)
        self.fileMenu.insertSeparator()
        self.fileExitAction.addTo(self.fileMenu)
        self.fileMenu.insertSeparator()
        self.MenuBar.insertItem(QString(""),self.fileMenu,1)

        self.Edit = QPopupMenu(self)
        self.editedit_atomsAction.addTo(self.Edit)
        self.editperiodic_expandAction.addTo(self.Edit)
        self.editgenerate_linkatomsAction.addTo(self.Edit)
        self.editHCS_link_layersAction.addTo(self.Edit)
        self.editBCTC_link_LayersAction.addTo(self.Edit)
        self.MenuBar.insertItem(QString(""),self.Edit,2)

        self.Statistics = QPopupMenu(self)
        self.statisticsCoordinationsAction.addTo(self.Statistics)
        self.statisticsChargesAction.addTo(self.Statistics)
        self.Statistics.insertSeparator()
        self.statisticssave_BCT_coefficientsAction.addTo(self.Statistics)
        self.statisticsget_charge_constraintsAction.addTo(self.Statistics)
        self.Statistics.insertSeparator()
        self.popupMenu_20 = QPopupMenu(self)
        self.Statistics.setAccel(QString.null,self.Statistics.insertItem(self.statisticsgraphsAction.iconSet(),self.__tr("graphs"),self.popupMenu_20))
        self.statisticsgraphselement_charge_histogramsAction.addTo(self.popupMenu_20)
        self.statisticsgraphsradial_distribution_functionAction.addTo(self.popupMenu_20)
        self.Statistics.insertSeparator()
        self.statisticsbond_listAction.addTo(self.Statistics)
        self.Statistics.insertSeparator()
        self.statisticsSave_StatisticsAction.addTo(self.Statistics)
        self.MenuBar.insertItem(QString(""),self.Statistics,3)

        self.helpMenu = QPopupMenu(self)
        self.helpContentsAction.addTo(self.helpMenu)
        self.helpAboutAction.addTo(self.helpMenu)
        self.MenuBar.insertItem(QString(""),self.helpMenu,4)


        self.languageChange()

        self.resize(QSize(530,440).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.fileOpenAction,SIGNAL("activated()"),self.fileOpen)
        self.connect(self.fileSaveAction,SIGNAL("activated()"),self.fileSave)
        self.connect(self.fileSaveAsAction,SIGNAL("activated()"),self.fileSaveAs)
        self.connect(self.fileExitAction,SIGNAL("activated()"),self.close)
        self.connect(self.helpContentsAction,SIGNAL("activated()"),self.helpContents)
        self.connect(self.helpAboutAction,SIGNAL("activated()"),self.helpAbout)
        self.connect(self.statisticsCoordinationsAction,SIGNAL("activated()"),self.DisplayCoordinations)
        self.connect(self.statisticsSave_StatisticsAction,SIGNAL("activated()"),self.SaveStatistics)
        self.connect(self.statisticsbond_listAction,SIGNAL("activated()"),self.display_bondlist)
        self.connect(self.editedit_atomsAction,SIGNAL("activated()"),self.editAtoms)
        self.connect(self.editAtomsButton,SIGNAL("clicked()"),self.editAtoms)
        self.connect(self.periodicExpandButton,SIGNAL("clicked()"),self.periodicExpand)
        self.connect(self.editperiodic_expandAction,SIGNAL("activated()"),self.periodicExpand)
        self.connect(self.statisticsChargesAction,SIGNAL("activated()"),self.DisplayCharges)
        self.connect(self.editgenerate_linkatomsAction,SIGNAL("activated()"),self.generate_linkatoms)
        self.connect(self.statisticsgraphselement_charge_histogramsAction,SIGNAL("activated()"),self.elementChargesHistograms)
        self.connect(self.statisticsgraphsradial_distribution_functionAction,SIGNAL("activated()"),self.rdf)
        self.connect(self.helpContentsAction,SIGNAL("activated()"),self.helpContents)
        self.connect(self.editBCTC_link_LayersAction,SIGNAL("activated()"),self.generate_BCTCLinkAtoms)
        self.connect(self.editHCS_link_layersAction,SIGNAL("activated()"),self.generate_hcsLinkAtoms)
        self.connect(self.statisticssave_BCT_coefficientsAction,SIGNAL("activated()"),self.save_BCTC_coefficients)
        self.connect(self.statisticsget_charge_constraintsAction,SIGNAL("activated()"),self.display_charge_constraints)

        self.setTabOrder(self.textBrowser1,self.editAtomsButton)
        self.setTabOrder(self.editAtomsButton,self.periodicExpandButton)
        self.setTabOrder(self.periodicExpandButton,self.whateverButton)
        self.setTabOrder(self.whateverButton,self.VMDButton)


    def languageChange(self):
        self.setCaption(self.__tr("Geometry Statistics and Editing"))
        QToolTip.add(self,self.__tr("Geometry Statitsics and Editing"))
        QWhatsThis.add(self,self.__tr("Tool to calculate gemeotry statistics and edit geometries"))
        QToolTip.add(self.frame3,self.__tr("Geometry operations quick access"))
        QWhatsThis.add(self.frame3,self.__tr("Geometry operations quick access"))
        self.periodicExpandButton.setText(self.__tr("Periodic Expand"))
        QToolTip.add(self.periodicExpandButton,self.__tr("periodically expand supercell"))
        QWhatsThis.add(self.periodicExpandButton,self.__tr("Opens a dialog to periodically expand the supercell"))
        self.editAtomsButton.setText(self.__tr("Edit Atoms"))
        QToolTip.add(self.editAtomsButton,self.__tr("edit atoms"))
        QWhatsThis.add(self.editAtomsButton,self.__tr("Opens a dialog to edit atom properties and add/remove geometry layers"))
        self.VMDButton.setText(self.__tr("VMD"))
        QToolTip.add(self.VMDButton,self.__tr("view geometry in VMD"))
        QWhatsThis.add(self.VMDButton,self.__tr("saves the current geometry and loads it into vmd"))
        self.whateverButton.setText(self.__tr("nothing yet"))
        QToolTip.add(self.whateverButton,self.__tr("nothing"))
        QWhatsThis.add(self.whateverButton,self.__tr("nothing"))
        self.textBrowser1.setText(self.__tr("geostats.py -- tool to analyze and manupulate molecular geometries"))
        QToolTip.add(self.textBrowser1,self.__tr("Textual geometry Statistics"))
        QWhatsThis.add(self.textBrowser1,self.__tr("Displays textual statistics information on the loaded geometry. Can be saved in .html format via Statistics->save statistics"))
        self.fileOpenAction.setText(self.__tr("Open"))
        self.fileOpenAction.setMenuText(self.__tr("&Open..."))
        self.fileOpenAction.setAccel(self.__tr("Ctrl+O"))
        self.fileSaveAction.setText(self.__tr("Save"))
        self.fileSaveAction.setMenuText(self.__tr("&Save"))
        self.fileSaveAction.setAccel(self.__tr("Ctrl+S"))
        self.fileSaveAsAction.setText(self.__tr("Save As"))
        self.fileSaveAsAction.setMenuText(self.__tr("Save &As..."))
        self.fileSaveAsAction.setAccel(QString.null)
        self.fileExitAction.setText(self.__tr("Exit"))
        self.fileExitAction.setMenuText(self.__tr("E&xit"))
        self.fileExitAction.setAccel(QString.null)
        self.helpContentsAction.setText(self.__tr("Contents"))
        self.helpContentsAction.setMenuText(self.__tr("&Contents..."))
        self.helpContentsAction.setAccel(QString.null)
        self.helpAboutAction.setText(self.__tr("About"))
        self.helpAboutAction.setMenuText(self.__tr("&About"))
        self.helpAboutAction.setAccel(QString.null)
        self.statisticsCoordinationsAction.setText(self.__tr("coordinations"))
        self.statisticsCoordinationsAction.setMenuText(self.__tr("coordinations"))
        self.statisticsSave_StatisticsAction.setText(self.__tr("save Statistics"))
        self.statisticsSave_StatisticsAction.setMenuText(self.__tr("save Statistics"))
        self.statisticsbond_listAction.setText(self.__tr("bond list"))
        self.statisticsbond_listAction.setMenuText(self.__tr("bond list"))
        self.statisticsallAction.setText(self.__tr("all"))
        self.statisticsallAction.setMenuText(self.__tr("all"))
        self.editedit_atomsAction.setText(self.__tr("edit atoms..."))
        self.editedit_atomsAction.setMenuText(self.__tr("edit atoms..."))
        self.editperiodic_expandAction.setText(self.__tr("periodic expand..."))
        self.editperiodic_expandAction.setMenuText(self.__tr("periodic expand..."))
        self.fileView_in_VMDAction.setText(self.__tr("View in VMD"))
        self.fileView_in_VMDAction.setMenuText(self.__tr("View in VMD"))
        self.statisticsChargesAction.setText(self.__tr("charges"))
        self.editgenerate_linkatomsAction.setText(self.__tr("simple-link layers..."))
        self.editgenerate_linkatomsAction.setMenuText(self.__tr("simple-link layers..."))
        self.statisticsgraphsAction.setText(self.__tr("graphs"))
        self.statisticsgraphsAction.setMenuText(self.__tr("graphs"))
        self.statisticsgraphselement_charge_histogramsAction.setText(self.__tr("element charge histograms..."))
        self.statisticsgraphselement_charge_histogramsAction.setMenuText(self.__tr("element charge histograms..."))
        self.statisticsgraphsradial_distribution_functionAction.setText(self.__tr("radial distribution function..."))
        self.statisticsgraphsradial_distribution_functionAction.setMenuText(self.__tr("radial distribution function..."))
        self.statisticssave_BCT_coefficientsAction.setText(self.__tr("save BCT coefficients..."))
        self.statisticssave_BCT_coefficientsAction.setMenuText(self.__tr("save BCT coefficients..."))
        self.statisticsget_charge_constraintsAction.setText(self.__tr("get charge constraints..."))
        self.statisticsget_charge_constraintsAction.setMenuText(self.__tr("get charge constraints..."))
        self.editBCTC_link_LayersAction.setText(self.__tr("BCTC-link layers..."))
        self.editBCTC_link_LayersAction.setMenuText(self.__tr("BCTC-link layers..."))
        self.editBCTC_link_LayersAction.setToolTip(self.__tr("Create QM/MM linked geometry using BCTC embedding"))
        self.editHCS_link_layersAction.setText(self.__tr("HCS-link layers..."))
        self.editHCS_link_layersAction.setMenuText(self.__tr("HCS-link layers..."))
        self.editHCS_link_layersAction.setToolTip(self.__tr("Generate QM/MM embedded geometry using HCS"))
        if self.MenuBar.findItem(1):
            self.MenuBar.findItem(1).setText(self.__tr("&File"))
        if self.MenuBar.findItem(2):
            self.MenuBar.findItem(2).setText(self.__tr("Edit"))
        self.Statistics.changeItem(self.Statistics.idAt(6),self.__tr("graphs"))
        if self.MenuBar.findItem(3):
            self.MenuBar.findItem(3).setText(self.__tr("Statistics"))
        if self.MenuBar.findItem(4):
            self.MenuBar.findItem(4).setText(self.__tr("&Help"))


    def fileOpen(self):
        	fn=str(QFileDialog.getOpenFileName(".","Known formats (*.gen *.xyz *.fmg)", self))
        	if (len(fn)!=0):
        		self.fileread(fn)
        

    def fileSave(self):
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
        

    def fileSaveAs(self):
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
        

    def DisplayCoordinations(self):
        	self.textBrowser1.append(self.geo.rt_coordinations())
        	self.statisticsCoordinationsAction.setDisabled(1)
        

    def helpContents(self):
        	gsh=geostatshelp.GeostatsHelp()
        	gsh.show()
        

    def helpAbout(self):
        	QMessageBox.about(self,"About geostats","""<H1>geostats %s</H1><p>Calculate geometry statistics and transform file formats.</p>
        			  <p>Licensed under the Non-Profit Open Software License version 3.0</p>
        			  <p>see file LICENSE for details.</p>
        			  <p>Written by Jan M. Knaup <b><a href="mailto:Jan.Knaup@bccms.uni-bremen.de">Jan.Knaup@bccms.uni-bremen.de</a></b></p>"""%(comatsci.constants.VERSION))
        

    def SaveStatistics(self):
        	statsfile=open("geostats.htm","w")
        	statsfile.write(str(self.textBrowser1.text()))
        

    def fileread(self,a0):
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
        

    def display_bondlist(self):
        	self.textBrowser1.append(self.geo.rt_bondlist())
        	self.statisticsbond_listAction.setDisabled(1)
        

    def editAtoms(self):
        	GeoEditor=EditGeometry.Edit_Atoms()
        	GeoEditor.setGeometry(self.geo)
        	GeoEditor.exec_loop()
        	self.geo=GeoEditor.originalgeo
        	if self.geo.layerbyname("PCHR")!=None:
        		self.editgenerate_linkatomsAction.setEnabled(True)
        

    def periodicExpand(self):
        	pedialog=periodicExpand.periodicExpand()
        	pedialog.setVectors(self.geo.Lattice)
        	if pedialog.exec_loop()==QDialog.Accepted:
        		self.geo.periodicexpand(pedialog.getExpandSpecifier())
        

    def DisplayCharges(self):
        	self.textBrowser1.append(self.geo.rt_charges())
        	self.statisticsChargesAction.setDisabled(1)
        

    def generate_linkatoms(self):
        	linkDialog=linkLayersSLA.layerLinkSLA()
        	linkDialog.setGeometry(self.geo)
        	if linkDialog.exec_loop()==QDialog.Accepted:
        		self.geo=linkDialog.embeddedGeometry
        

    def generate_hcsLinkAtoms(self):
        	linkDialog=linkLayersHCS.layerLinkHCS()
        	linkDialog.setGeometry(self.geo)
        	if linkDialog.exec_loop()==QDialog.Accepted:
        		self.geo=linkDialog.embeddedGeometry
        

    def elementChargesHistograms(self):
        	elcharhist=elementChargeHistograms.elementChargesHistograms()
        	elcharhist.setGeometry(self.geo)
        	elcharhist.exec_loop()
        

    def rdf(self):
        	rdfwiz=rdfWizard.rdfWizard()
        	rdfwiz.setGeometry(self.geo)
        	rdfwiz.exec_loop()
        

    def generate_BCTCLinkAtoms(self):
        	linkDialog=linkLayersBCTC.layerLinkBCTC()
        	linkDialog.setGeometry(self.geo)
        	if linkDialog.exec_loop()==QDialog.Accepted:
        		self.geo=linkDialog.embeddedGeometry
        

    def save_BCTC_coefficients(self):
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
        

    def display_charge_constraints(self):
        	# we have a wizard for this, just for convenience
        	constraintsDialog=chargeConstraintsWizard.atomChargeConstraintsWizard()
        	# pass the geomerty object to the wizard
        	constraintsDialog.setGeometry(self.geo)
        	# The wizard does everything by itself, so just run it
        	constraintsDialog.exec_loop()
        

    def __tr(self,s,c = None):
        return qApp.translate("MainWindow",s,c)
