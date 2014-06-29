# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'geostatsHelp.ui'
#
# Created: Thu Sep 27 11:28:57 2012
#      by: PyQt4 UI code generator 4.7.2
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_GeostatsHelp(object):
    def setupUi(self, GeostatsHelp):
        GeostatsHelp.setObjectName("GeostatsHelp")
        GeostatsHelp.resize(486, 480)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/images/images/GSLogoHelp.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        GeostatsHelp.setWindowIcon(icon)
        self.hboxlayout = QtGui.QHBoxLayout(GeostatsHelp)
        self.hboxlayout.setObjectName("hboxlayout")
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.textBrowser = QtGui.QTextBrowser(GeostatsHelp)
        self.textBrowser.setObjectName("textBrowser")
        self.verticalLayout.addWidget(self.textBrowser)
        self.hboxlayout.addLayout(self.verticalLayout)

        self.retranslateUi(GeostatsHelp)
        QtCore.QMetaObject.connectSlotsByName(GeostatsHelp)

    def retranslateUi(self, GeostatsHelp):
        GeostatsHelp.setWindowTitle(QtGui.QApplication.translate("GeostatsHelp", "Geostats Help", None, QtGui.QApplication.UnicodeUTF8))
        GeostatsHelp.setWindowIconText(QtGui.QApplication.translate("GeostatsHelp", "Help on Geostats", None, QtGui.QApplication.UnicodeUTF8))
        GeostatsHelp.setToolTip(QtGui.QApplication.translate("GeostatsHelp", "Help on geostats.py.py", None, QtGui.QApplication.UnicodeUTF8))
        GeostatsHelp.setWhatsThis(QtGui.QApplication.translate("GeostatsHelp", "Help on geostats.py.py", None, QtGui.QApplication.UnicodeUTF8))
        self.textBrowser.setDocumentTitle(QtGui.QApplication.translate("GeostatsHelp", "Geostats Help", None, QtGui.QApplication.UnicodeUTF8))
        self.textBrowser.setHtml(QtGui.QApplication.translate("GeostatsHelp", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><title>Geostats Help</title><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:18px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:xx-large; font-weight:600;\">Geostats Help</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><img src=\":/images/images/GSLogo.svg\" /></p></body></html>", None, QtGui.QApplication.UnicodeUTF8))

import gsresources_rc
