from PyQt5.uic import loadUiType
import PyQt5.QtCore as QtCore

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import pdb
Ui_MainWindow, QMainWindow = loadUiType('./mainWindow.ui')
Ui_SignalWindow, QSignalWindow = loadUiType('./signal.ui')
class Main(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super(Main, self).__init__()
        self.signals = {}
        self.materials = {}
        self.matRefs = []
        self.orientations = {}
        self.orientRefs = []
        
        self.setupUi(self)
        self.signalButton.clicked.connect(self.setupSignal)
        self.matPropButton.clicked.connect(self.setupMatProp)
        
    def addGeometryPlot(self,fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
    
    def rmGeometryPlot(self):
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        
    def setupSignal(self):
        
class SignalWindow(QSignalWindow,Ui_SignalWindow):
    def __init__(self):
        super(SignalWindow,self).__init__()
        self.setupUi(self)
        self.inputNodeTable.installEventFilter(self)
        
    def eventFilter(self,source, event):    
        if event.type() == QtCore.QEvent.KeyPress:
            if event.key() == QtCore.Qt.Key_V and event.modifiers() == QtCore.Qt.ControlModifier:
                clip = QtWidgets.QApplication.clipboard()
                mime = clip.mimeData()
                data = mime.data('application/x-qt-windows-mime;value="Csv"')
                data = str(data.data())
                data = data.split('\r\n')
                source.setRowCount(source.rowCount()+len(data)-1)
                
                for i,row in enumerate(data[:-1]):
                    rowData = np.fromstring(row,sep=',')
                    source.setItem(i,0,QtWidgets.QTableWidgetItem(str(rowData[0])))
                    source.setItem(i,1,QtWidgets.QTableWidgetItem(str(rowData[1])))
                return True
###############Following is not working properly yet                
            elif event.key() == QtCore.Qt.Key_Right and source.currentColumn()<source.columnCount():
                print [source.currentColumn(), source.currentRow()]
                source.setCurrentCell(source.currentRow(),source.currentColumn()+1)
                print [source.currentColumn(), source.currentRow()]
            elif event.key() == QtCore.Qt.Key_Left and source.currentColumn()>0:
                print [source.currentColumn(), source.currentRow()]
                source.setCurrentCell(source.currentRow(),source.currentColumn()-1)
                print [source.currentColumn(), source.currentRow()]
            elif event.key() == QtCore.Qt.Key_Up and source.currentRow()>0:
                print [source.currentColumn(), source.currentRow()]
                source.setCurrentCell(source.currentRow()-1,source.currentColumn())
                print [source.currentColumn(), source.currentRow()]
            elif event.key() == QtCore.Qt.Key_Down and  source.currentRow()<source.rowCount():
                print [source.currentColumn(), source.currentRow()]
                source.setCurrentCell(source.currentRow()+1,source.currentColumn())
                print [source.currentColumn(), source.currentRow()]
        return False
        
        

            
if __name__ == '__main__':
    import sys
    from PyQt5 import QtWidgets
    import numpy as np
    fig1 = Figure()
    ax1f1 = fig1.add_subplot(111)
    ax1f1.plot(np.random.rand(5))
    app = QtWidgets.QApplication(sys.argv)
    main = Main()
    main.addGeometryPlot(fig1)
    main.show()
    sys.exit(app.exec_())