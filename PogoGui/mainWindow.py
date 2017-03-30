# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainWindow.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(955, 647)
        MainWindow.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.ImageWindow = QtWidgets.QWidget(self.centralwidget)
        self.ImageWindow.setGeometry(QtCore.QRect(9, 9, 501, 451))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ImageWindow.sizePolicy().hasHeightForWidth())
        self.ImageWindow.setSizePolicy(sizePolicy)
        self.ImageWindow.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.ImageWindow.setObjectName("ImageWindow")
        self.mplvl = QtWidgets.QVBoxLayout(self.ImageWindow)
        self.mplvl.setContentsMargins(0, 0, 0, 0)
        self.mplvl.setObjectName("mplvl")
        self.signalButton = QtWidgets.QPushButton(self.centralwidget)
        self.signalButton.setGeometry(QtCore.QRect(600, 130, 75, 23))
        self.signalButton.setObjectName("signalButton")
        self.seamButton = QtWidgets.QPushButton(self.centralwidget)
        self.seamButton.setGeometry(QtCore.QRect(600, 160, 75, 23))
        self.seamButton.setObjectName("seamButton")
        self.matPropButton = QtWidgets.QPushButton(self.centralwidget)
        self.matPropButton.setGeometry(QtCore.QRect(544, 190, 131, 23))
        self.matPropButton.setObjectName("matPropButton")
        self.orientButton = QtWidgets.QPushButton(self.centralwidget)
        self.orientButton.setGeometry(QtCore.QRect(544, 220, 131, 23))
        self.orientButton.setObjectName("orientButton")
        self.BCButton = QtWidgets.QPushButton(self.centralwidget)
        self.BCButton.setGeometry(QtCore.QRect(544, 250, 131, 23))
        self.BCButton.setObjectName("BCButton")
        self.ElementSize = QtWidgets.QTextEdit(self.centralwidget)
        self.ElementSize.setGeometry(QtCore.QRect(570, 300, 104, 21))
        self.ElementSize.setPlaceholderText("")
        self.ElementSize.setObjectName("ElementSize")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(540, 280, 131, 20))
        self.label.setObjectName("label")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 955, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.signalButton.setText(_translate("MainWindow", "Add signal"))
        self.seamButton.setText(_translate("MainWindow", "Add seam"))
        self.matPropButton.setText(_translate("MainWindow", "Add material properties"))
        self.orientButton.setText(_translate("MainWindow", "Add orientation"))
        self.BCButton.setText(_translate("MainWindow", "Add boundary condition"))
        self.label.setText(_translate("MainWindow", "Approximate Element Size"))

