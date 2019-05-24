# -*- coding: utf-8 -*-
import sys
import math
import numpy as np
import math_part
import matplotlib.pyplot as plt
from SplForm import *
# Импортируем наш интерфейс из файла
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QWidget, QTableWidget, QTableWidgetItem, QHeaderView
from numpy import float64
from PyQt5 import QtWidgets, QtGui, QtCore
from matplotlib.figure import Figure
from tabwidg import *
class MyWin(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, *args, **kwargs):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.pushButton.clicked.connect(self.MyFunction)
        self.comboBox.currentTextChanged.connect(self.Changed)


    def MyFunction(self):

        def f_plot(*args, **kwargs):
            xlist = []
            ylist = []
            for i, arg in enumerate(args):
                if(i % 2 == 0):
                    xlist.append(arg)
                else:
                    ylist.append(arg) 
    
            colors = kwargs.pop('colors', 'k')
            linewidth = kwargs.pop('linewidth', 1.)
    
            fig = plt.figure()
            ax = fig.add_subplot(111)
            i = 0
            for x, y, color in zip(xlist, ylist, colors):
                i += 1
                ax.plot(x, y, color=color, linewidth=linewidth, label=str(i))
    
            ax.grid(True)
            ax.legend()
            plt.show()
            #save('ex_1_6_1', fmt='pdf')
            #save('ex_1_6_1', fmt='png')

        if (self.comboBox.currentText() == 'ТЕСТОВАЯ'):
            flag=0
        elif (self.comboBox.currentText() == 'ОСНОВНАЯ 1'):
            flag=1
        elif (self.comboBox.currentText() == 'ОСНОВНАЯ 2'):
            flag=2
        elif (self.comboBox.currentText() == 'ОСНОВНАЯ 2.1'):
            flag=3
        if(self.checkBox.isChecked()):
            natSpl=1
        else: natSpl=0

        n = int(self.textEdit.toPlainText())
        cnt=int(self.textEdit_2.toPlainText())
        N = n*cnt
        A = float64(self.textEdit_3.toPlainText())
        B = float64(self.textEdit_4.toPlainText())
        h=float64((B-A)/n)
        hs=float64(h/cnt)
        x=np.zeros(n+1, float64)
        y=np.zeros(n+1, float64)
        xN=np.zeros(N+1, float64)
        yN=np.zeros(N+1, float64)
        dyN=np.zeros(N+1, float64)
        ddyN=np.zeros(N+1, float64)

        x,y=math_part.Grid(x,y, A,B,h,n, flag)
        xN,yN=math_part.Grid(xN,yN,A,B,hs,N, flag)
        for i in range(N+1):
            dyN[i]=math_part.df(xN[i],flag)
            ddyN[i]=math_part.ddf(xN[i],flag)
        
        a,b,c,d = math_part.calculate_coefficients(x,y,n, N, A, B,h, natSpl, flag)

        spline=np.zeros(N+1, float64)
        splineFirstDev=np.zeros(N+1, float64)
        splineSecDev=np.zeros(N+1, float64)
        p=1
        for i in range(N+1):
            if(x[p]<xN[i]):
                p+=1
            spline[i]=math_part.s(x[p],xN[i],a[p],b[p],c[p],d[p])
            splineFirstDev[i]=math_part.ds(x[p],xN[i],a[p],b[p],c[p],d[p])
            splineSecDev[i]=math_part.dds(x[p],xN[i],a[p],b[p],c[p],d[p])

        self.Table(x, y,n,a,b,c,d)
        self.Table2(xN, yN,spline,splineFirstDev,dyN,N,a,b,c,d)
        self.Table3(xN, ddyN, splineSecDev,dyN,N)

        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        colors = [ 'blue','red']
        if(self.checkBox_2.isChecked()):
            f_plot(x, y,xN, spline,colors=colors,linewidth=2.)
            f_plot(xN, dyN,xN, splineFirstDev,colors=colors,linewidth=2.)
            f_plot(xN, ddyN,xN, splineSecDev,colors=colors,linewidth=2.)
        max1, max2,max3,ix1,ix2,ix3=math_part.Spravka(self,yN, spline, dyN, splineFirstDev, ddyN, splineSecDev, N)
        #SPRAVKA
        x1=xN[ix1]
        x2=xN[ix2]
        x3=xN[ix3]
        self.Label(n, N, max1, max2, max3, x1,x2,x3)
        # plt.grid(True)
        # plt.show()
        # self.secWin = second_window(self)        
        #self.secWin.show()
    def Label(self,n, N,max1, max2, max3, x1,x2,x3):
        self.textEdit_5.setText(
            "Сетка сплайна:  n = " +str(N)+'\n'+
            "Основная сетка: N = " + str(n) +'\n'+
            "max|F(xj)-S(xj)| = "  +str(max1) + "    при x = "+str(x1)+'\n'
            "max|F'(xj)-S'(xj)| = "  +str(max2) + "    при x = "+str(x2)+'\n'
            "max|F''(xj)-S''(xj)| = "  +str(max3) + "    при x = "+str(x3)+'\n'
        )

    def Table(self, x, y, n,a,b,c,d):
        self.tableWidget.setRowCount(n)
        for i in range(0, n):
            self.tableWidget.setItem(i, 0, QtWidgets.QTableWidgetItem(str(i+1)))
            self.tableWidget.setItem(i, 1, QtWidgets.QTableWidgetItem(str(x[i])))
            self.tableWidget.setItem(i, 2, QtWidgets.QTableWidgetItem(str(x[i+1])))
            self.tableWidget.setItem(i, 3, QtWidgets.QTableWidgetItem(str(a[i+1])))
            self.tableWidget.setItem(i, 4, QtWidgets.QTableWidgetItem(str(b[i+1])))
            self.tableWidget.setItem(i, 5, QTableWidgetItem(str(c[i+1])))
            self.tableWidget.setItem(i, 6, QTableWidgetItem(str(d[i+1])))
        self.tableWidget.resizeColumnsToContents()

    def Table2(self, xN, yN,s,ds,dyN,N,a,b,c,d):
        self.tableWidget_2.setRowCount(N+1)
        for i in range(0, N+1):
            self.tableWidget_2.setItem(i, 0, QtWidgets.QTableWidgetItem(str(i)))
            self.tableWidget_2.setItem(i, 1, QtWidgets.QTableWidgetItem(str(xN[i])))
            self.tableWidget_2.setItem(i, 2, QtWidgets.QTableWidgetItem(str(yN[i])))
            self.tableWidget_2.setItem(i, 3, QtWidgets.QTableWidgetItem(str(s[i])))
            self.tableWidget_2.setItem(i, 4, QtWidgets.QTableWidgetItem(str(abs(yN[i]-s[i]))))
            self.tableWidget_2.setItem(i, 5, QTableWidgetItem(str(dyN[i])))
            self.tableWidget_2.setItem(i, 6, QTableWidgetItem(str(ds[i])))
            self.tableWidget_2.setItem(i, 7, QTableWidgetItem(str(abs(dyN[i]-ds[i]))))
        self.tableWidget_2.resizeColumnsToContents()

    def Table3(self, xN, ddyN, dds,dyN,N):
        self.tableWidget_3.setRowCount(N+1)
        for i in range(0, N+1):
            self.tableWidget_3.setItem(i, 0, QtWidgets.QTableWidgetItem(str(i)))
            self.tableWidget_3.setItem(i, 1, QtWidgets.QTableWidgetItem(str(xN[i])))
            self.tableWidget_3.setItem(i, 2, QtWidgets.QTableWidgetItem(str(ddyN[i])))
            self.tableWidget_3.setItem(i, 3, QtWidgets.QTableWidgetItem(str(dds[i])))
            self.tableWidget_3.setItem(i, 4, QtWidgets.QTableWidgetItem(str(abs(ddyN[i]-dds[i]))))
        self.tableWidget_3.resizeColumnsToContents()

    def Changed(self):
        if (self.comboBox.currentText() == 'ТЕСТОВАЯ'):
            self.textEdit_3.setText(str(-1))
            self.textEdit_4.setText(str(1))
        else:
            self.textEdit_3.setText(str(2))
            self.textEdit_4.setText(str(4))

class second_window(QMainWindow, Ui_MainWindow_tab):
    def __init__(self, parent=None, *args, **kwargs):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)

if __name__=="__main__":
    app = QtWidgets.QApplication(sys.argv)
    myapp = MyWin()
    myapp.show()
    try:
        sys.exit(app.exec_())
    except SystemExit:
        pass

# # -*- coding: utf-8 -*-
# import sys
# #import math
# import math_part
# # Импортируем наш интерфейс из файла

# from SplForm import *
# from PyQt5.QtWidgets import QApplication, QMainWindow
# from MyMplCanc import MtMplCanv
# from numpy import float64
# from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
# from PyQt5 import QtWidgets, QtGui, QtCore
# from MyMplCanc import MtMplCanv
# from MyMplCanc import MtMplCanv2
# from matplotlib.figure import Figure
# from tabwidg import *

# class MyWin(QMainWindow, Ui_MainWindow):

#     def __init__(self, parent=None):
#         QMainWindow.__init__(self, parent)
#         self.setupUi(self)
#         self.secWin = None
#         self.figure = Figure()
#         self.figure2 = Figure()

#         # добавление шаблона размещения на виджет
#         self.companovka_for_mpl = QtWidgets.QVBoxLayout(self.widget)
#         # получение объекта класса холста с нашим рисунком

#         self.canvas = MtMplCanv(self.figure)
#         # Размещение экземпляра класса холста в шаблоне размещения
        
#         self.companovka_for_mpl.addWidget(self.canvas)
#         # получение объекта класса панели управления холста
        
#         self.toolbar = NavigationToolbar(self.canvas, self)
#         # Размещение экземпляра класса панели управления в шаблоне размещения
        
#         self.companovka_for_mpl.addWidget(self.toolbar)
#         # Здесь прописываем событие нажатия на кнопку
        
#         self.pushButton.clicked.connect(self.MyFunction)
#     def MyFunction(self):
#         if str(self.comboBox.currentText()) == 'ТЕСТОВАЯ':
#             # self.textEdit_3.toPlainText('-1')
#             # self.textEdit_4.toPlainText('1')
#             A=-1
#             B=1
#         else: 
#             A=2
#             B=4
        
#         # if str(self.comboBox.currentText()) == 'ОСНОВНАЯ 1':
#         #     self.textBrowser.setText("f1")
        
#         # if str(self.comboBox.currentText()) == 'ОСНОВНАЯ 2':
#         #     self.textBrowser.setText("f2")

#         # if str(self.comboBox.currentText()) == 'ОСНОВНАЯ 2.1':
#         #     self.textBrowser.setText("f2.1")

#         n    = int(self.textEdit.toPlainText())
#         nDiv = int(self.textEdit_2.toPlainText())
        
#         # A = float64(self.textEdit_3.toPlainText())
#         # B = float64(self.textEdit_4.toPlainText())
        
#         self.secWin = second_window(self)
#         math_part.mathpart.building(self, n, nDiv, A, B, self.secWin)
#         self.secWin.show()

# class second_window(QMainWindow, Ui_MainWindow_tab):
#     def __init__(self, parent=None, *args, **kwargs):
#         QMainWindow.__init__(self, parent)
#         self.setupUi(self)


# if __name__=="__main__":
#     app = QtWidgets.QApplication(sys.argv)
#     myapp = MyWin()
#     myapp.show()
#     try:
#         sys.exit(app.exec_())
#     except SystemExit:
#         pass