# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'signal.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(272, 388)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setGeometry(QtCore.QRect(30, 220, 161, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.AmplitudeSelector = QtWidgets.QPushButton(Dialog)
        self.AmplitudeSelector.setGeometry(QtCore.QRect(50, 190, 121, 23))
        self.AmplitudeSelector.setObjectName("AmplitudeSelector")
        self.inputNodeTable = QtWidgets.QTableWidget(Dialog)
        self.inputNodeTable.setGeometry(QtCore.QRect(20, 20, 211, 141))
        self.inputNodeTable.setObjectName("inputNodeTable")
        self.inputNodeTable.setColumnCount(2)
        self.inputNodeTable.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.inputNodeTable.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.inputNodeTable.setHorizontalHeaderItem(1, item)

        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.AmplitudeSelector.setText(_translate("Dialog", "Choose amplitude file"))
        item = self.inputNodeTable.horizontalHeaderItem(0)
        item.setText(_translate("Dialog", "New Column"))
        item = self.inputNodeTable.horizontalHeaderItem(1)
        item.setText(_translate("Dialog", "New Column"))

