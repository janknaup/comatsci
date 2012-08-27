# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dftbChargeConstraintDialog.ui'
#
# Created: Mon Aug 27 17:48:19 2012
#      by: PyQt4 UI code generator 4.7.2
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(563, 480)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/images/images/GSLogo.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Dialog.setWindowIcon(icon)
        self.verticalLayoutWidget = QtGui.QWidget(Dialog)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(0, 0, 561, 481))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_2 = QtGui.QLabel(self.verticalLayoutWidget)
        self.label_2.setWordWrap(True)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.serialListPlainTextEdit = QtGui.QPlainTextEdit(self.verticalLayoutWidget)
        self.serialListPlainTextEdit.setTabChangesFocus(True)
        self.serialListPlainTextEdit.setObjectName("serialListPlainTextEdit")
        self.verticalLayout.addWidget(self.serialListPlainTextEdit)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.preafactorLabel = QtGui.QLabel(self.verticalLayoutWidget)
        self.preafactorLabel.setObjectName("preafactorLabel")
        self.horizontalLayout.addWidget(self.preafactorLabel)
        self.prefactorSpinBox = QtGui.QDoubleSpinBox(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.prefactorSpinBox.sizePolicy().hasHeightForWidth())
        self.prefactorSpinBox.setSizePolicy(sizePolicy)
        self.prefactorSpinBox.setMinimum(0.0)
        self.prefactorSpinBox.setMaximum(9999999.0)
        self.prefactorSpinBox.setSingleStep(0.05)
        self.prefactorSpinBox.setProperty("value", 1.0)
        self.prefactorSpinBox.setObjectName("prefactorSpinBox")
        self.horizontalLayout.addWidget(self.prefactorSpinBox)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.constraintOutputTextBrowser = QtGui.QTextBrowser(self.verticalLayoutWidget)
        self.constraintOutputTextBrowser.setTabChangesFocus(True)
        self.constraintOutputTextBrowser.setTextInteractionFlags(QtCore.Qt.TextSelectableByKeyboard|QtCore.Qt.TextSelectableByMouse)
        self.constraintOutputTextBrowser.setOpenLinks(False)
        self.constraintOutputTextBrowser.setObjectName("constraintOutputTextBrowser")
        self.verticalLayout.addWidget(self.constraintOutputTextBrowser)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.constraintButton = QtGui.QPushButton(self.verticalLayoutWidget)
        self.constraintButton.setObjectName("constraintButton")
        self.horizontalLayout_2.addWidget(self.constraintButton)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.closeButton = QtGui.QPushButton(self.verticalLayoutWidget)
        self.closeButton.setObjectName("closeButton")
        self.horizontalLayout_2.addWidget(self.closeButton)
        self.verticalLayout.addLayout(self.horizontalLayout_2)

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.closeButton, QtCore.SIGNAL("released()"), Dialog.accept)
        QtCore.QMetaObject.connectSlotsByName(Dialog)
        Dialog.setTabOrder(self.serialListPlainTextEdit, self.prefactorSpinBox)
        Dialog.setTabOrder(self.prefactorSpinBox, self.constraintButton)
        Dialog.setTabOrder(self.constraintButton, self.constraintOutputTextBrowser)
        Dialog.setTabOrder(self.constraintOutputTextBrowser, self.closeButton)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "DFTB+ charge constraints", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("Dialog", "Input serial numbers (counting from 1) of atoms to apply constraints on, separated by white space, in the text editor below. \n"
"Click \"constraints\" to generate constraints string.", None, QtGui.QApplication.UnicodeUTF8))
        self.preafactorLabel.setText(QtGui.QApplication.translate("Dialog", "charge constraints prefactor", None, QtGui.QApplication.UnicodeUTF8))
        self.constraintButton.setText(QtGui.QApplication.translate("Dialog", "constraints", None, QtGui.QApplication.UnicodeUTF8))
        self.closeButton.setText(QtGui.QApplication.translate("Dialog", "close", None, QtGui.QApplication.UnicodeUTF8))

import gsresources_rc
