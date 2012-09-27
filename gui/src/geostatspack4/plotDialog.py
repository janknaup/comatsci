# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'plotDialog.ui'
#
# Created: Thu Sep 27 11:28:58 2012
#      by: PyQt4 UI code generator 4.7.2
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(631, 606)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/images/images/GSLogo.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Dialog.setWindowIcon(icon)
        self.frame_2 = QtGui.QFrame(Dialog)
        self.frame_2.setGeometry(QtCore.QRect(10, 510, 621, 41))
        self.frame_2.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtGui.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.horizontalLayoutWidget = QtGui.QWidget(self.frame_2)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(0, 0, 621, 41))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtGui.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.addOutputCheckBox = QtGui.QCheckBox(self.horizontalLayoutWidget)
        self.addOutputCheckBox.setChecked(True)
        self.addOutputCheckBox.setObjectName("addOutputCheckBox")
        self.horizontalLayout.addWidget(self.addOutputCheckBox)
        self.saveDataCheckBox = QtGui.QCheckBox(self.horizontalLayoutWidget)
        self.saveDataCheckBox.setObjectName("saveDataCheckBox")
        self.horizontalLayout.addWidget(self.saveDataCheckBox)
        self.label = QtGui.QLabel(self.horizontalLayoutWidget)
        self.label.setEnabled(False)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.dataFileLineEdit = QtGui.QLineEdit(self.horizontalLayoutWidget)
        self.dataFileLineEdit.setEnabled(False)
        self.dataFileLineEdit.setObjectName("dataFileLineEdit")
        self.horizontalLayout.addWidget(self.dataFileLineEdit)
        self.outputGroupBox = QtGui.QGroupBox(Dialog)
        self.outputGroupBox.setGeometry(QtCore.QRect(10, 0, 611, 341))
        self.outputGroupBox.setMinimumSize(QtCore.QSize(611, 341))
        self.outputGroupBox.setMaximumSize(QtCore.QSize(611, 341))
        self.outputGroupBox.setObjectName("outputGroupBox")
        self.parametersGroupBox = QtGui.QGroupBox(Dialog)
        self.parametersGroupBox.setGeometry(QtCore.QRect(10, 350, 611, 141))
        self.parametersGroupBox.setMinimumSize(QtCore.QSize(611, 141))
        self.parametersGroupBox.setMaximumSize(QtCore.QSize(611, 141))
        self.parametersGroupBox.setObjectName("parametersGroupBox")
        self.horizontalLayoutWidget_2 = QtGui.QWidget(Dialog)
        self.horizontalLayoutWidget_2.setGeometry(QtCore.QRect(10, 560, 621, 41))
        self.horizontalLayoutWidget_2.setObjectName("horizontalLayoutWidget_2")
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.horizontalLayoutWidget_2)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.plotButton = QtGui.QPushButton(self.horizontalLayoutWidget_2)
        self.plotButton.setObjectName("plotButton")
        self.horizontalLayout_2.addWidget(self.plotButton)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.saveButton = QtGui.QPushButton(self.horizontalLayoutWidget_2)
        self.saveButton.setObjectName("saveButton")
        self.horizontalLayout_2.addWidget(self.saveButton)
        self.cancelButton = QtGui.QPushButton(self.horizontalLayoutWidget_2)
        self.cancelButton.setObjectName("cancelButton")
        self.horizontalLayout_2.addWidget(self.cancelButton)

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.saveDataCheckBox, QtCore.SIGNAL("toggled(bool)"), self.dataFileLineEdit.setEnabled)
        QtCore.QObject.connect(self.saveDataCheckBox, QtCore.SIGNAL("toggled(bool)"), self.label.setEnabled)
        QtCore.QObject.connect(self.cancelButton, QtCore.SIGNAL("released()"), Dialog.reject)
        QtCore.QObject.connect(self.saveButton, QtCore.SIGNAL("released()"), Dialog.accept)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        self.addOutputCheckBox.setText(QtGui.QApplication.translate("Dialog", "add to output", None, QtGui.QApplication.UnicodeUTF8))
        self.saveDataCheckBox.setText(QtGui.QApplication.translate("Dialog", "save data", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("Dialog", "data filename suffix", None, QtGui.QApplication.UnicodeUTF8))
        self.outputGroupBox.setTitle(QtGui.QApplication.translate("Dialog", "Output", None, QtGui.QApplication.UnicodeUTF8))
        self.parametersGroupBox.setTitle(QtGui.QApplication.translate("Dialog", "Parameters", None, QtGui.QApplication.UnicodeUTF8))
        self.plotButton.setText(QtGui.QApplication.translate("Dialog", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.saveButton.setText(QtGui.QApplication.translate("Dialog", "Save", None, QtGui.QApplication.UnicodeUTF8))
        self.cancelButton.setText(QtGui.QApplication.translate("Dialog", "Cancel", None, QtGui.QApplication.UnicodeUTF8))

import gsresources_rc
