# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'geostatsMain.ui'
#
# Created: Thu Sep 27 11:28:58 2012
#      by: PyQt4 UI code generator 4.7.2
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(632, 547)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/images/images/GSLogo.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.textBrowser = QtGui.QTextBrowser(self.centralwidget)
        self.textBrowser.setProperty("cursor", QtCore.Qt.ArrowCursor)
        self.textBrowser.setObjectName("textBrowser")
        self.verticalLayout.addWidget(self.textBrowser)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 632, 23))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuEdit = QtGui.QMenu(self.menubar)
        self.menuEdit.setObjectName("menuEdit")
        self.menuStatistics = QtGui.QMenu(self.menubar)
        self.menuStatistics.setObjectName("menuStatistics")
        self.menuSummaries = QtGui.QMenu(self.menuStatistics)
        self.menuSummaries.setEnabled(False)
        self.menuSummaries.setObjectName("menuSummaries")
        self.menuGraphs = QtGui.QMenu(self.menuStatistics)
        self.menuGraphs.setEnabled(False)
        self.menuGraphs.setObjectName("menuGraphs")
        self.menuHelp = QtGui.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionOpen = QtGui.QAction(MainWindow)
        self.actionOpen.setObjectName("actionOpen")
        self.actionSave_geometry = QtGui.QAction(MainWindow)
        self.actionSave_geometry.setEnabled(False)
        self.actionSave_geometry.setObjectName("actionSave_geometry")
        self.actionSave_geometry_as = QtGui.QAction(MainWindow)
        self.actionSave_geometry_as.setEnabled(False)
        self.actionSave_geometry_as.setObjectName("actionSave_geometry_as")
        self.actionSave_output = QtGui.QAction(MainWindow)
        self.actionSave_output.setObjectName("actionSave_output")
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName("actionExit")
        self.actionAbout = QtGui.QAction(MainWindow)
        self.actionAbout.setObjectName("actionAbout")
        self.actionHelp = QtGui.QAction(MainWindow)
        self.actionHelp.setObjectName("actionHelp")
        self.actionEdit_Atoms = QtGui.QAction(MainWindow)
        self.actionEdit_Atoms.setEnabled(False)
        self.actionEdit_Atoms.setObjectName("actionEdit_Atoms")
        self.actionPeriodic_Expand = QtGui.QAction(MainWindow)
        self.actionPeriodic_Expand.setEnabled(False)
        self.actionPeriodic_Expand.setObjectName("actionPeriodic_Expand")
        self.actionSimple_link_Layers = QtGui.QAction(MainWindow)
        self.actionSimple_link_Layers.setEnabled(False)
        self.actionSimple_link_Layers.setObjectName("actionSimple_link_Layers")
        self.actionBCTC_link_Layers = QtGui.QAction(MainWindow)
        self.actionBCTC_link_Layers.setEnabled(False)
        self.actionBCTC_link_Layers.setObjectName("actionBCTC_link_Layers")
        self.actionCoordinations = QtGui.QAction(MainWindow)
        self.actionCoordinations.setObjectName("actionCoordinations")
        self.actionCharges = QtGui.QAction(MainWindow)
        self.actionCharges.setObjectName("actionCharges")
        self.actionWrite_BCTC_Coefficients = QtGui.QAction(MainWindow)
        self.actionWrite_BCTC_Coefficients.setEnabled(False)
        self.actionWrite_BCTC_Coefficients.setObjectName("actionWrite_BCTC_Coefficients")
        self.actionDFTB_Charge_Constraints = QtGui.QAction(MainWindow)
        self.actionDFTB_Charge_Constraints.setEnabled(False)
        self.actionDFTB_Charge_Constraints.setObjectName("actionDFTB_Charge_Constraints")
        self.actionRadial_Distribution_Functions = QtGui.QAction(MainWindow)
        self.actionRadial_Distribution_Functions.setObjectName("actionRadial_Distribution_Functions")
        self.actionAngle_Distribution_Histogram = QtGui.QAction(MainWindow)
        self.actionAngle_Distribution_Histogram.setObjectName("actionAngle_Distribution_Histogram")
        self.actionBond_Length_Histogram = QtGui.QAction(MainWindow)
        self.actionBond_Length_Histogram.setObjectName("actionBond_Length_Histogram")
        self.actionElement_Charge_Histograms = QtGui.QAction(MainWindow)
        self.actionElement_Charge_Histograms.setObjectName("actionElement_Charge_Histograms")
        self.actionPrint_Output = QtGui.QAction(MainWindow)
        self.actionPrint_Output.setObjectName("actionPrint_Output")
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionSave_geometry)
        self.menuFile.addAction(self.actionSave_geometry_as)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionPrint_Output)
        self.menuFile.addAction(self.actionSave_output)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menuEdit.addAction(self.actionEdit_Atoms)
        self.menuEdit.addAction(self.actionPeriodic_Expand)
        self.menuEdit.addSeparator()
        self.menuEdit.addAction(self.actionSimple_link_Layers)
        self.menuEdit.addAction(self.actionBCTC_link_Layers)
        self.menuSummaries.addAction(self.actionCoordinations)
        self.menuSummaries.addAction(self.actionCharges)
        self.menuGraphs.addAction(self.actionRadial_Distribution_Functions)
        self.menuGraphs.addAction(self.actionAngle_Distribution_Histogram)
        self.menuGraphs.addAction(self.actionBond_Length_Histogram)
        self.menuGraphs.addSeparator()
        self.menuGraphs.addAction(self.actionElement_Charge_Histograms)
        self.menuStatistics.addAction(self.menuSummaries.menuAction())
        self.menuStatistics.addAction(self.menuGraphs.menuAction())
        self.menuStatistics.addSeparator()
        self.menuStatistics.addAction(self.actionWrite_BCTC_Coefficients)
        self.menuStatistics.addAction(self.actionDFTB_Charge_Constraints)
        self.menuHelp.addAction(self.actionAbout)
        self.menuHelp.addAction(self.actionHelp)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())
        self.menubar.addAction(self.menuStatistics.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QObject.connect(self.actionExit, QtCore.SIGNAL("activated()"), MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "geostats.py", None, QtGui.QApplication.UnicodeUTF8))
        self.textBrowser.setDocumentTitle(QtGui.QApplication.translate("MainWindow", "geostats.py", None, QtGui.QApplication.UnicodeUTF8))
        self.textBrowser.setHtml(QtGui.QApplication.translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><title>geostats.py</title><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p align=\"center\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><img src=\":/images/images/GSLogo.svg\" /></p>\n"
"<p align=\"center\" style=\" margin-top:18px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:xx-large; font-weight:600;\">geostats.py</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">tool to analyze and manupulate molecular geometries for computational materials science</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">part of the comatsci computational materials science toolkit</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setTitle(QtGui.QApplication.translate("MainWindow", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.menuEdit.setTitle(QtGui.QApplication.translate("MainWindow", "Edit", None, QtGui.QApplication.UnicodeUTF8))
        self.menuStatistics.setTitle(QtGui.QApplication.translate("MainWindow", "Statistics", None, QtGui.QApplication.UnicodeUTF8))
        self.menuSummaries.setTitle(QtGui.QApplication.translate("MainWindow", "Summaries", None, QtGui.QApplication.UnicodeUTF8))
        self.menuGraphs.setTitle(QtGui.QApplication.translate("MainWindow", "Graphs", None, QtGui.QApplication.UnicodeUTF8))
        self.menuHelp.setTitle(QtGui.QApplication.translate("MainWindow", "Help", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen.setText(QtGui.QApplication.translate("MainWindow", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+O", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave_geometry.setText(QtGui.QApplication.translate("MainWindow", "Save Geometry", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave_geometry.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+S", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave_geometry_as.setText(QtGui.QApplication.translate("MainWindow", "Save Geometry As", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave_geometry_as.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Shift+S", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave_output.setText(QtGui.QApplication.translate("MainWindow", "Save Output", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave_output.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+W", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setText(QtGui.QApplication.translate("MainWindow", "Exit", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+Q", None, QtGui.QApplication.UnicodeUTF8))
        self.actionAbout.setText(QtGui.QApplication.translate("MainWindow", "About", None, QtGui.QApplication.UnicodeUTF8))
        self.actionHelp.setText(QtGui.QApplication.translate("MainWindow", "Help", None, QtGui.QApplication.UnicodeUTF8))
        self.actionHelp.setShortcut(QtGui.QApplication.translate("MainWindow", "F1", None, QtGui.QApplication.UnicodeUTF8))
        self.actionEdit_Atoms.setText(QtGui.QApplication.translate("MainWindow", "Edit Atoms", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPeriodic_Expand.setText(QtGui.QApplication.translate("MainWindow", "Periodic Expand", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSimple_link_Layers.setText(QtGui.QApplication.translate("MainWindow", "Simple-link Layers", None, QtGui.QApplication.UnicodeUTF8))
        self.actionBCTC_link_Layers.setText(QtGui.QApplication.translate("MainWindow", "BCTC-link Layers", None, QtGui.QApplication.UnicodeUTF8))
        self.actionCoordinations.setText(QtGui.QApplication.translate("MainWindow", "Coordinations", None, QtGui.QApplication.UnicodeUTF8))
        self.actionCharges.setText(QtGui.QApplication.translate("MainWindow", "Charges", None, QtGui.QApplication.UnicodeUTF8))
        self.actionWrite_BCTC_Coefficients.setText(QtGui.QApplication.translate("MainWindow", "Write BCTC Coefficients", None, QtGui.QApplication.UnicodeUTF8))
        self.actionDFTB_Charge_Constraints.setText(QtGui.QApplication.translate("MainWindow", "DFTB Charge Constraints", None, QtGui.QApplication.UnicodeUTF8))
        self.actionRadial_Distribution_Functions.setText(QtGui.QApplication.translate("MainWindow", "Radial Distribution Functions", None, QtGui.QApplication.UnicodeUTF8))
        self.actionAngle_Distribution_Histogram.setText(QtGui.QApplication.translate("MainWindow", "Angle Distribution Histogram", None, QtGui.QApplication.UnicodeUTF8))
        self.actionBond_Length_Histogram.setText(QtGui.QApplication.translate("MainWindow", "Bond Length Histogram", None, QtGui.QApplication.UnicodeUTF8))
        self.actionElement_Charge_Histograms.setText(QtGui.QApplication.translate("MainWindow", "Element Charge Histograms", None, QtGui.QApplication.UnicodeUTF8))
        self.actionElement_Charge_Histograms.setToolTip(QtGui.QApplication.translate("MainWindow", "Element Charge Histograms", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPrint_Output.setText(QtGui.QApplication.translate("MainWindow", "Print Output", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPrint_Output.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+P", None, QtGui.QApplication.UnicodeUTF8))

import gsresources_rc
