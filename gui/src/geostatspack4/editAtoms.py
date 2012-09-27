# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'editAtoms.ui'
#
# Created: Thu Sep 27 11:28:59 2012
#      by: PyQt4 UI code generator 4.7.2
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(680, 710)
        Dialog.setMinimumSize(QtCore.QSize(680, 710))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/images/images/GSLogo.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Dialog.setWindowIcon(icon)
        Dialog.setSizeGripEnabled(False)
        Dialog.setModal(False)
        self.verticalLayout_2 = QtGui.QVBoxLayout(Dialog)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.groupBox = QtGui.QGroupBox(Dialog)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout = QtGui.QVBoxLayout(self.groupBox)
        self.verticalLayout.setObjectName("verticalLayout")
        self.atomsTableView = QtGui.QTableView(self.groupBox)
        self.atomsTableView.setObjectName("atomsTableView")
        self.verticalLayout.addWidget(self.atomsTableView)
        self.frame_2 = QtGui.QFrame(self.groupBox)
        self.frame_2.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtGui.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.frame_2)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.addAtomButton = QtGui.QPushButton(self.frame_2)
        self.addAtomButton.setObjectName("addAtomButton")
        self.horizontalLayout_2.addWidget(self.addAtomButton)
        self.delAtomButton = QtGui.QPushButton(self.frame_2)
        self.delAtomButton.setObjectName("delAtomButton")
        self.horizontalLayout_2.addWidget(self.delAtomButton)
        self.verticalLayout.addWidget(self.frame_2)
        self.verticalLayout_2.addWidget(self.groupBox)
        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Reset)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout_2.addWidget(self.buttonBox)
        self.verticalLayout_2.setStretch(0, 10)
        self.verticalLayout_2.setStretch(1, 1)

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("accepted()"), Dialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("rejected()"), Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Edit Atoms", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox.setTitle(QtGui.QApplication.translate("Dialog", "Atoms", None, QtGui.QApplication.UnicodeUTF8))
        self.addAtomButton.setText(QtGui.QApplication.translate("Dialog", "Add Atom", None, QtGui.QApplication.UnicodeUTF8))
        self.delAtomButton.setText(QtGui.QApplication.translate("Dialog", "Delete Atom", None, QtGui.QApplication.UnicodeUTF8))

import gsresources_rc
