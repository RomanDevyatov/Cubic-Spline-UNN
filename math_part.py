import math
import numpy as np
from numpy import float64
from matplotlib.figure import Figure
from SplForm import Ui_MainWindow
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5 import QtWidgets, QtGui, QtCore
import matplotlib.pyplot as plt


def f(x, flag):
    if(flag==0):
        if (x<=0):
            return math.pow(x,3) + 3.0 * math.pow(x,2)
        elif (x > 0):
            return -math.pow(x,3) + 3.0 * math.pow(x,2)
    if(flag==1):
        return math.log(x + 1) / x
    if(flag==2):
        return math.log(x + 1) / x + math.cos(10*x)
    if(flag==3):
        return math.log(x + 1) / x + math.cos(100*x)

def df(x, flag):
    if(flag==0):    
        if (x<=0):
            return 3.0*math.pow(x,2) + 6.0 * x
        elif (x > 0):
            return -3.0*math.pow(x,2) + 6.0 * x
    if(flag==1):
        return 1.0 / (x * (1 + x)) - math.log(x + 1) / (x * x)
    if(flag==2):
        return 1.0 / (x * (1 + x)) - math.log(x + 1) / (x * x) - 10.0 * math.sin(10 * x)
    if(flag==3):
        return 1.0 / (x * (1 + x)) - math.log(x + 1) / (x * x) - 100.0 * math.sin(100 * x)

def ddf(x, flag):
    if(flag==0):
        if (x<=0):
            return 6.0*x + 6.0
        elif (x > 0):
            return -6.0*x + 6.0
    if(flag==1):
        return -1.0 / (x * (1 + x) * (1 + x)) - 2.0 / (math.pow(x,2) * (1 + x)) + 2.0 * math.log(x + 1) / math.pow(x,3)
    if(flag==2):
        return -1.0 / (x * (1 + x) * (1 + x)) - 2.0 / (math.pow(x,2) * (1 + x)) + 2.0 * math.log(x + 1) / math.pow(x,3) - 100.0 * math.cos(10 * x)
    if(flag==3):
        return -1.0 / (x * (1 + x) * (1 + x)) - 2.0 / (math.pow(x,2) * (1 + x)) + 2.0 * math.log(x + 1) / math.pow(x,3) - 10000.0 * math.cos(100 * x)

def s(x,xCurr,a,b,c,d):
    return a+b*(xCurr-x)+c/2.0*math.pow((xCurr-x), 2)+d/6.0*math.pow(xCurr-x,3)

def ds(x,xCurr,a,b,c,d):
    return b+c*(xCurr-x)+d/2.0*math.pow(xCurr-x,2)

def dds(x,xCurr,a,b,c,d):
    return c+d*(xCurr-x)

def TDMASolve(self,mu1,mu2,a,b,c,n,h):
    alpha=np.zeros(n+1)
    beta=np.zeros(n+1)

    alpha[1]=0
    beta[1]=mu1
    for i in range(1,n,1):
        alpha[i + 1] = h / (-4 * h - alpha[i] * h)
        beta[i + 1] = ((-6.0 / h) * (a[i + 1] - 2 * a[i] + a[i - 1]) + beta[i] * h) / (-4 * h - alpha[i] * h)
    c[n]=mu2
    #for i in range(n,0,-1):
    i=n
    while i>0:
        c[i-1]=alpha[i]*c[i]+beta[i]
        i=i-1
    return c

def TDMASolvelog(a, c, b, d):
    n = len(d)-1
    alpha = np.zeros(n)
    beta = np.zeros(n)
    x = np.zeros(n+1)
    alpha[0] = 0 # alpha[0] = kapa1, kapa1 в этой задаче всегда 0
    beta[0] = d[0] # beta[0] = mu[1]
    for i in range(0, n-1):
        alpha[i+1] = b[i]/(c[i]-alpha[i]*a[i])
        beta[i+1] = (d[i+1] + beta[i] * a[i])/(c[i] - alpha[i] * a[i])

    x[-1] = d[-1]
    for i in range (n-1, 0, -1):
        x[i] = alpha[i] * x[i+1] + beta[i]
    x[0] = alpha[0] * x[1] + beta[0]
    return x

def calculate_coefficients(x,y,n, N, A, B,h, natSplFlag, flag):
    a=np.zeros(n+1, float64)
    b=np.zeros(n+1, float64)
    c=np.zeros(n+1, float64)
    d=np.zeros(n+1, float64) 

    for i in range(n+1):        
        a[i] = f(A+i*h, flag)
    mu1=0
    mu2=0
    if(natSplFlag==0):
        mu1 = ddf(A,flag)
        mu2 = ddf(B, flag)
    
    #c=TDMASolve(mu1,mu2,a,b,c,n,h)
    #############################
    # alpha=np.zeros(n+1)
    # beta=np.zeros(n+1)

    # alpha[1]=0
    # beta[1]=mu1
    # for i in range(1,n,1):
    #     alpha[i + 1] = h / (-4 * h - alpha[i] * h)
    #     beta[i + 1] = ((-6.0 / h) * (a[i + 1] - 2 * a[i] + a[i - 1]) + beta[i] * h) / (-4 * h - alpha[i] * h)
    # c[n]=mu2
    # for i in range(n,0,-1):
    #     c[i-1]=alpha[i]*c[i]+beta[i]
    #################################
    C = np.zeros(n-1, float64)
    A = np.zeros(n-1, float64)
    B = np.zeros(n-1, float64)
    # ????????????????
    d= np.zeros(n+1, float64)
    for i in range(1, n):  # последний индекс - n-1
            A[i-1] = 1
            C[i-1] = -4
            B[i-1] = 1
    d[0]=mu1
        
    for i in range(1, n):
            d[i] = -6*(y[i-1] - 2*y[i] + y[i+1])/(h**2)

    nl = len(d)-1
    alpha = np.zeros(nl)
    beta = np.zeros(nl)
    x = np.zeros(nl+1)
    alpha[0] = 0 # alpha[0] = kapa1, kapa1 в этой задаче всегда 0
    beta[0] = d[0] # beta[0] = mu[1]
    for i in range(0, nl-1):
        alpha[i+1] = b[i]/(c[i]-alpha[i]*a[i])
        beta[i+1] = (d[i+1] + beta[i] * a[i])/(c[i] - alpha[i] * a[i])

    x[-1] = d[-1]
    for i in range (nl-1, 0, -1):
        x[i] = alpha[i] * x[i+1] + beta[i]
    x[0] = alpha[0] * x[1] + beta[0]
    c=x

    #???????????????
    for i in range(1,n+1):
        d[i] = (c[i] - c[i-1]) / h
        b[i] = (a[i] - a[i - 1]) / h + h*c[i]/3.0 + h*c[i-1]/ 6.0
            
    return a,b,c,d


def Grid(x,y, A,B,h,n,flag):
    for i in range(n):
        x[i]=A+i*h
        y[i]=f(x[i],flag)
    x[n]=B
    y[n]=f(x[n],flag)

    return x, y


def Spravka(self,yN, s, dyN, ds, ddyN, dds, N):
    maxYS=0
    maxDYDS=0
    maxDDYDDS=0
    i1=0
    i2=0
    i3=0
    for i in range(N+1):
        if(maxYS<(abs(yN[i]-s[i]))):
            maxYS=(abs(yN[i]-s[i]))
            i1=i
        if(maxDYDS<(abs(dyN[i]-ds[i]))):
            maxDYDS=(abs(dyN[i]-ds[i]))
            i2=i
        if(maxDDYDDS<(abs(ddyN[i]-dds[i]))):
            maxDDYDDS=(abs(ddyN[i]-dds[i]))
            i3=i
    return maxYS, maxDYDS, maxDDYDDS, i1, i2, i3





# class mathpart(Ui_MainWindow):

    

    # def calculate_spline(self, n):
    #     def function_F(x):
    #         if self.comboBox.currentText() == "Тестовая":
    #             return x**3 + 3*x**2 if x <= 0 else -x**3 + 3*x**2
    #         elif self.comboBox.currentText() == "log(x+1) / x":
    #             return math.log1p(x) / x
    #         elif self.comboBox.currentText() == "log(x+1) / x + cos(10x)":
    #             return math.log1p(x) / x + math.cos(10*x)
    #         elif self.comboBox.currentText() == "log(x+1) / x + cos(100x)":
    #             return math.log1p(x) / x + math.cos(100*x)

    #     def derivative_F(x):
    #         if self.comboBox.currentText() == "Тестовая":
    #             return 3*x**2 + 6*x if x <= 0 else -3*x**2 + 6*x
    #         elif self.comboBox.currentText() == "log(x+1) / x":
    #             return 1.0 / (x * (1 + x)) - math.log1p(x) / (x**2)
    #         elif self.comboBox.currentText() == "log(x+1) / x + cos(10x)":
    #             return 1.0 / (x * (1 + x)) - math.log1p(x) / (x **2) - 10 * math.sin(10 * x)
    #         elif self.comboBox.currentText() == "log(x+1) / x + cos(100x)":
    #             return 1.0 / (x * (1 + x)) - math.log1p(x) / (x **2) - 100 * math.sin(100 * x)

    #     def second_der_F(x):
    #         if self.comboBox.currentText() == "Тестовая":
    #             return 6*x+6 if x <= 0 else -6*x + 6
    #         elif self.comboBox.currentText() == "log(x+1) / x":
    #             return -1.0 / (x * (1 + x)**2) - 2.0 / (x**2 * (1 + x)) + 2.0 * math.log1p(x) / x**3
    #         elif self.comboBox.currentText() == "log(x+1) / x + cos(10x)":
    #             return -1.0 / (x * (x + 1)**2) - 2.0 / (x**2 * (1 + x)) + 2.0 * math.log1p(x) / x**3 - 100 * math.cos(10 * x)
    #         elif self.comboBox.currentText() == "log(x+1) / x + cos(100x)":
    #             return -1.0 / (x * (x + 1)**2) - 2.0 / (x**2 * (1 + x)) + 2.0 * math.log1p(x) / x**3 - 10000 * math.cos(100 * x)
    #     def S(x, xAdd, a, b , c, d):                        
    #         return a+b*(xAdd - x)+c/2.0*(xAdd - x)**2+d/6.0 * (xAdd - x)**3

    #     def dS(x, xAdd, a, b , c, d):
    #         return b+c*(xAdd-x)+d/2.0*(xAdd-x)**2

    #     def ddS(x, xAdd, a, b , c, d):
    #         return c+d*(xAdd-x)

        
    #     left = -1 if self.comboBox.currentText() == "Тестовая" else 2
    #     right = 1 if self.comboBox.currentText() == "Тестовая" else 4
    #     h = (right-left)/n 
    #     c = np.zeros(n+1, float64)
    #     g = np.zeros(n+1, float64)
    #     C = np.zeros(n-1, float64)
    #     A = np.zeros(n-1, float64)
    #     B = np.zeros(n-1, float64)
    #     for i in range(1, n):  # последний индекс - n-1
    #         A[i-1] = 1
    #         C[i-1] = -4
    #         B[i-1] = 1
    #     if self.comboBox_2.currentText() == "Совпадение 2 произв.":
    #         g[0] = second_der_F(left)
    #     elif self.comboBox_2.currentText() == "ЕГУ":
    #         g[0] = 0
    #     f = np.zeros(n+1, float64)
    #     for i in range(0, n+1):
    #         f[i] = function_F(left + i*h)
    #     for i in range(1, n):
    #         g[i] = -6*(f[i-1] - 2*f[i] + f[i+1])/(h**2)
    #     if self.comboBox_2.currentText() == "Совпадение 2 произв.":
    #         g[0] = second_der_F(right)
    #     elif self.comboBox_2.currentText() == "ЕГУ":
    #         g[0] = 0

    #     c = TDMASolve(A, C, B, g)
    #     b = np.zeros(n, float64)
    #     d = np.zeros(n, float64)
    #     a = np.zeros(n, float64)
        
    #     for i in range(0, n):
    #         d[i] = (c[i+1] - c[i])/h
    #         a[i] = f[i+1]
    #         b[i] = (f[i+1] - f[i])/h + (2*c[i+1] + c[i]) * h / 6
    #     c = c[1:]
    #     self.tableWidget.setRowCount(n)
    #     for i in range(0, n):
    #         self.tableWidget.setItem(i, 0, QtWidgets.QTableWidgetItem(str(i)))
    #         self.tableWidget.setItem(
    #             i, 1, QtWidgets.QTableWidgetItem(str(left + i*h)))
    #         self.tableWidget.setItem(
    #             i, 2, QtWidgets.QTableWidgetItem(str(left + (i+1)*h)))
    #         self.tableWidget.setItem(
    #             i, 3, QtWidgets.QTableWidgetItem(str(a[i])))
    #         self.tableWidget.setItem(
    #             i, 4, QtWidgets.QTableWidgetItem(str(b[i])))
    #         self.tableWidget.setItem(
    #             i, 5, QtWidgets.QTableWidgetItem(str(c[i])))
    #         self.tableWidget.setItem(
    #             i, 6, QtWidgets.QTableWidgetItem(str(d[i])))

    #     plt.subplot(111)
    #     x = left
    #     y = np.zeros(4*n, float64)
    #     for i in range(0, n):
    #         x_i = left+(i+1)*h
    #         y[4*i] = (a[i] + b[i] * (x-x_i) + c[i]*((x-x_i)**2)/2 + d[i]*((x-x_i)**3)/6)
    #         x += h/4
    #         y[4*i+1] =  (a[i] + b[i] * (x-x_i) + c[i]*((x-x_i)**2)/2 + d[i]*((x-x_i)**3)/6)
    #         x += h/4
    #         y[4*i+2] =  (a[i] + b[i] * (x-x_i) + c[i]*((x-x_i)**2)/2 + d[i]*((x-x_i)**3)/6)
    #         x += h/4
    #         y[4*i+3] =  (a[i] + b[i] * (x-x_i) + c[i]*((x-x_i)**2)/2 + d[i]*((x-x_i)**3)/6)
    #         x += h/4
    #     plt.plot(y)

    #     x = left
    #     y_der = []
    #     for i in range(0, n):
    #         x_i = left+(i+1)*h
    #         y_der.append(b[i] + c[i]*(x-x_i) + d[i]/2*(x-x_i)**2)
    #         x += h/4
    #         y_der.append(b[i] + c[i]*(x-x_i) + d[i]/2*(x-x_i)**2)
    #         x += h/4
    #         y_der.append(b[i] + c[i]*(x-x_i) + d[i]/2*(x-x_i)**2)
    #         x += h/4
    #         y_der.append(b[i] + c[i]*(x-x_i) + d[i]/2*(x-x_i)**2)
    #         x += h/4
    #     plt.plot(y_der)
    #     x = left
    #     y_sec_der = []
    #     for i in range(0, n):
    #         x_i = left +(i+1)*h
    #         y_sec_der.append(c[i] + d[i]*(x-x_i))
    #         x += h/4
    #         y_sec_der.append(c[i] + d[i]*(x-x_i))
    #         x += h/4
    #         y_sec_der.append(c[i] + d[i]*(x-x_i))
    #         x += h/4
    #         y_sec_der.append(c[i] + d[i]*(x-x_i))
    #         x += h/4
    #     plt.plot(y_sec_der)

    #     N = 4*n
    #     h = (right-left)/(N)

    #     f = np.zeros(N, float64)
    #     for i in range(N):
    #         f[i] = function_F(left + i*h)
    #     f_der = np.zeros(N, float64)
    #     for i in range(N):
    #         f_der[i] = derivative_F(left + i*h)
    #     f_sec_der = np.zeros(N, float64)
    #     for i in range(N):
    #         f_sec_der[i] = second_der_F(left + i*h)
    #     plt.plot(f)
    #     plt.plot(f_der)
    #     plt.plot(f_sec_der)
    #     plt.legend(("Сплайн", "1 производная", "2 производная", "Функция", "1 производная", "2 производная"))
    #     plt.show()


    #     #spravka######################################################
    #     t=0
    #     maxFS=0
    #     maxdFdS=0
    #     maxddFddS=0
    #     fs_x=0
    #     dfds_x=0
    #     ddfdds_x=0
    #     xAdd=left
    #     while t<=N:
    #         k=1
    #         if(xAdd+t*h/4 > xAdd+k*h):
    #             k+=1
    #         if( maxFS>abs(function_F(xAdd+t*h/4)-S(xAdd+k*h, xAdd+t*h/4, a[k],b[k],c[k],d[k])) ):
    #             maxFS=abs(function_F(xAdd+t*h/4)-S(xAdd+k*h, xAdd+t*h/4, a[k],b[k],c[k],d[k]))
    #             fs_x=xAdd+t*h/4

    #         if( maxdFdS>abs(derivative_F(xAdd+t*h/4)-dS(xAdd+k*h, xAdd+t*h/4, a[k],b[k],c[k],d[k])) ):
    #             maxdFdS=abs(derivative_F(xAdd+t*h/4)-dS(xAdd+k*h, xAdd+t*h/4, a[k],b[k],c[k],d[k]))
    #             dfds_x=xAdd+t*h/4

    #         if( maxddFddS>abs(second_der_F(xAdd+t*h/4)-ddS(xAdd+k*h, xAdd+t*h/4, a[k],b[k],c[k],d[k])) ):
    #             maxddFddS=abs(second_der_F(xAdd+t*h/4)-ddS(xAdd+k*h, xAdd+t*h/4, a[k],b[k],c[k],d[k]))
    #             ddfdds_x=xAdd+t*h/4
    #         t+=1
        
    #     # secwin.textBrowser.setText("max|F(xj)− S(xj)| = " + str(maxFS) + " при x = " + str(fs_x))
    #     # secwin.textBrowser_2.setText("max|F'(xj)− S'(xj)| = " + str(maxdFdS)+ " при x = " + str(dfds_x))
    #     # secwin.textBrowser_3.setText("max|F''(xj)− S''(xj)| = " + str(maxddFddS)+ " при x = " + str(ddfdds_x))
    #     ########################################################################################





    #     y_der = np.array(y_der)
    #     y_sec_der = np.array(y_sec_der)

    #     self.label.setText(QtCore.QCoreApplication.translate(
    #         "MainWindow",
    #         "Справка \nСетка сплайна: n = «" + str(n) + "» \n" +
    #         "Контрольная сетка: N = «" + str(N) + "» \n" +
    #         "Погрешность сплайна на контрольной сетке \n" +
    #         "max F(x) - S(x) = " + str(max(abs(y-f))) + " при x = " + str(-1 + np.argmax(abs(y-f))*h) + "\n" +
    #         "Погрешность производной на контрольной сетке \n" +
    #         "max F'(x) - S'(x) = " + str(max(abs(y_der-f_der))) + " при x = " + str(-1 + np.argmax(abs(y_der-f_der))*h) + "\n" +
    #         "Погрешность второй производной на контрольной сетке \n" +
    #         "max F''(x) - S''(x) = " + str(max(abs(y_sec_der-f_sec_der))) + " при x = " + str(-1 + np.argmax(abs(y_sec_der-f_sec_der))*h)))
    #     self.tableWidget_2.setRowCount(N)
    #     for i in range(0, N):
    #         self.tableWidget_2.setItem(i, 0, QtWidgets.QTableWidgetItem(str(i)))
    #         self.tableWidget_2.setItem(
    #             i, 1, QtWidgets.QTableWidgetItem(str(left + i*h)))
    #         self.tableWidget_2.setItem(
    #             i, 2, QtWidgets.QTableWidgetItem(str(f[i])))
    #         self.tableWidget_2.setItem(
    #             i, 3, QtWidgets.QTableWidgetItem(str(y[i])))
    #         self.tableWidget_2.setItem(
    #             i, 4, QtWidgets.QTableWidgetItem(str(f[i] - y[i])))
    #         self.tableWidget_2.setItem(
    #             i, 5, QtWidgets.QTableWidgetItem(str(f_der[i])))
    #         self.tableWidget_2.setItem(
    #             i, 6, QtWidgets.QTableWidgetItem(str(y_der[i])))
    #         self.tableWidget_2.setItem(
    #             i, 7, QtWidgets.QTableWidgetItem(str(f_der[i] - y_der[i])))
        
    #     plt.subplot(111)
    #     plt.plot(abs(y-f))
    #     plt.plot(abs(y_der-f_der))
    #     plt.plot(abs(y_sec_der-f_sec_der))
    #     plt.legend(("погрешность сплайна",
    #                 "погрешность 1 производной", "погрешность 2 производной"))
    #     plt.show()

# def TDMASolve(a, c, b, d):
#     n = len(d)-1
#     alpha = np.zeros(n)
#     beta = np.zeros(n)
#     x = np.zeros(n+1)
#     alpha[0] = 0 # alpha[0] = kapa1, kapa1 в этой задаче всегда 0
#     beta[0] = d[0] # beta[0] = mu[1]
#     for i in range(0, n-1):
#         alpha[i+1] = b[i]/(c[i]-alpha[i]*a[i])
#         beta[i+1] = (d[i+1] + beta[i] * a[i])/(c[i] - alpha[i] * a[i])

#     x[-1] = d[-1]
#     for i in range (n-1, 0, -1):
#         x[i] = alpha[i] * x[i+1] + beta[i]
#     x[0] = alpha[0] * x[1] + beta[0]
#     return x



# import math
# import pylab
# import numpy as np
# from numpy import float64
# from matplotlib import mlab
# from matplotlib.figure import Figure
# from SplForm import Ui_MainWindow
# from tab_widg import Ui_MainWindow_tab
# from PyQt5.QtWidgets import QApplication, QMainWindow
# from PyQt5 import QtWidgets, QtGui, QtCore
# from main import MyWin
# from main import second_window
# import matplotlib.pyplot as plt


# class mathpart(Ui_MainWindow):
#     def building(self, n, nDiv, A, B, secwin):

#         def f(x):
#             if(str(self.comboBox.currentText()) == 'ТЕСТОВАЯ'):
#                 if(x<=0):
#                     return x**3 + 3.0 * x**2
#                 elif (x > 0):
#                     return - x**3 + 3.0 * x**2
#             if(comboBoxBtn.currentText() == "ОСНОВНАЯ 1"):
#                 return Math.log1p(x) / x 
#             if(comboBoxBtn.currentText() == "ОСНОВНАЯ 2"):
#                 return Math.log1p(x) / x + Math.cos(10*x)
#             if(comboBoxBtn.currentText() == "ОСНОВНАЯ 2.1"):
#                 return Math.log1p(x) / x + Math.cos(100*x)
#         def df(x):
#             if(str(self.comboBox.currentText()) == 'ТЕСТОВАЯ'):
#                 if(x<=0):
#                     return 3.0 * x **2 + 6.0 * x
#                 elif (x > 0):
#                     return - 3.0 *x**2 + 6.0 * x
#             if(comboBoxBtn.currentText() == "ОСНОВНАЯ 1"):
#                 return 1.0 / (x * (1 + x)) - Math.log1p(x) / (x**2)
#             if(comboBoxBtn.currentText() == "ОСНОВНАЯ 2"):
#                 return  1.0 / (x * (1 + x)) - Math.log1p(x) / (x **2) - 10 * Math.sin(10 * x)
#             if(comboBoxBtn.currentText() == "ОСНОВНАЯ 2.1"):
#                 return 1.0 / (x * (1 + x)) - Math.log1p(x) / (x **2) - 100 * Math.sin(100 * x)
        
#         def ddf(x):
#             if(str(self.comboBox.currentText()) == 'ТЕСТОВАЯ'):
#                 if(x<=0):
#                     return 6 * x + 6
#                 elif (x > 0):
#                     return - 6 * x + 6
#             if(comboBoxBtn.currentText() == "ОСНОВНАЯ 1"):
#                 return -1.0 / (x * (1 + x)**2) - 2.0 / (x**2 * (1 + x)) + 2.0 * Math.log1p(x) / x**3
#             if(comboBoxBtn.currentText() == "ОСНОВНАЯ 2"):
#                 return -1.0 / (x * (x + 1)**2) - 2.0 / (x**2 * (1 + x)) + 2.0 * Math.log1p(x) / x**3 - 100 * Math.cos(10 * x)
#             if(comboBoxBtn.currentText() == "ОСНОВНАЯ 2.1"):
#                 return -1.0 / (x * (x + 1)**2) - 2.0 / (x**2 * (1 + x)) + 2.0 * Math.log1p(x) / x**3 - 10000 * Math.cos(100 * x)

#         def S(x, xAdd, a, b , c, d):                        
#             return a+b*(xAdd - x)+c/2.0*(xAdd - x)**2+d/6.0 * (xAdd - x)**3

#         def dS(x, xAdd, a, b , c, d):
#             return b+c*(xAdd-x)+d/2.0*(xAdd-x)**2

#         def ddS(x, xAdd, a, b , c, d):
#             return c+d*(xAdd-x)

#         def TDMASolve(a, c, b, d):
#             n = len(d)-1
#             alpha = np.zeros(n)
#             beta = np.zeros(n)
#             x = np.zeros(n+1)
#             alpha[0] = 0 # alpha[0] = kapa1, kapa1 в этой задаче всегда 0
#             beta[0] = d[0] # beta[0] = mu[1]
#             for i in range(0, n-1):
#                 alpha[i+1] = b[i]/(c[i]-alpha[i]*a[i])
#                 beta[i+1] = (d[i+1] + beta[i] * a[i])/(c[i] - alpha[i] * a[i])

#             x[-1] = d[-1]
#             for i in range (n-1, 0, -1):
#                 x[i] = alpha[i] * x[i+1] + beta[i]
#             x[0] = alpha[0] * x[1] + beta[0]
#             return x
        
#         N=n*nDiv

#         self.tableWidget.setRowCount(n)
#         self.tableWidget_2.setRowCount(N + 1)
#         self.tableWidget_3.setRowCount(N + 1)
        
#         x = np.zeros(n+1,float64)
#         y = np.zeros(n+1,float64)
#         xAdd = np.zeros(N+1)
#         yAdd = np.zeros(N+1)
#         h=(B-A)/n
#         sh=h/nDiv

#         #основная сетка
#         i = 0
#         while i < n:
#             x[i]= A+h*i
#             y[i]= f(x[i])
#             i += 1            
#         x[n]=B
#         y[n]=f(B)
                
        
#         #дополнительная сетка
#         i = 0
#         while i < N:
#             xAdd[i]=A+i*sh
#             yAdd[i]=f(xAdd[i])
#             i=i+1
#         xAdd[N]=B
#         yAdd[N]=f(B)
        
#         a=np.zeros(n+1)
#         b=np.zeros(n+1)
#         c=np.zeros(n+1)
#         d=np.zeros(n+1)

#         #коэффициент a
#         i = 0
#         while i <= n:
#             a[i]=f(x[i])
#             i=i+1
        
#         #функция F
#         F=np.zeros(n)
#         i = 2
#         while i <= n:            
#             F[i-1] = 6.0* ((a[i] - a[i-1]) / h - (a[i-1] - a[i-2]) / h)
#             i=i+1

#         # МЕТОД ПРОГОНКИ
#         g = np.zeros(n+1, float64)
#         C = np.zeros(n-1, float64)
#         A = np.zeros(n-1, float64)
#         B = np.zeros(n-1, float64)
#         left=A
#         for i in range(1, n):  # последний индекс - n-1
#             A[i-1] = 1
#             C[i-1] = -4
#             B[i-1] = 1
                    
#         if(self.checkBox.isChecked()):
#             g[0] = 0
#         else:
#             g[0] = ddf(A)
#         #f = np.zeros(n+1, float64)
#         # for i in range(0, n+1):
#         #     f[i] = function_F(left + i*h)
#         for i in range(1, n):
#             g[i] = -6*(y[i-1] - 2*y[i] + y[i+1])/(h**2)
#         if(self.checkBox.isChecked()):
#             g[0] = 0
#         else:
#             g[0] = ddf(B)

#         c = TDMASolve(A, C, B, g)
#         # alpha = np.zeros(n+2)
#         # beta = np.zeros(n+2)

#         # alpha[1]=0
#         # beta[1]=ddf(A)

#         # i=1
#         # while i <= n-1:
#         #     alpha[i + 1] = h / (-4 * h - alpha[i] * h)
#         #     beta[i + 1] = ((-6.0 / h) * (a[i + 1] - 2 * a[i] + a[i - 1]) + beta[i] * h) / (-4*h - alpha[i] * h)
#         #     i+=1
#         # c[n]=ddf(B)
#         # i=n
#         # while i >= 1:
#         #     c[i - 1] = alpha[i] * c[i] + beta[i]
#         #     i-=1

#         # if(self.checkBox.isChecked()):
#         #     c[0]=0
#         #     c[n]=0

#         i = 2
#         while i <= n:
#             print(str(F[i-1])+" = "+str(h*c[i-2]+2.0*(h+h)*c[i-1]+h*c[i]))
#             i=i+1

#         i = 1
#         while i <= n:
#             d[i] = (c[i] - c[i-1]) / h
#             b[i] = (a[i] - a[i - 1]) / h + h*c[i]/3.0 + h*c[i-1]/ 6.0
#             i+=1

#         #ТАБЛИЦЫ
#         i=1
#         while i<=n:
#             self.tableWidget.setItem(i, 0, QtWidgets.QTableWidgetItem(str(i)))
#             self.tableWidget.setItem(i, 1, QtWidgets.QTableWidgetItem(str(x[i-1])))
#             self.tableWidget.setItem(i, 2, QtWidgets.QTableWidgetItem(str(x[i])))
#             self.tableWidget.setItem(i, 3, QtWidgets.QTableWidgetItem(str(a[i])))
#             self.tableWidget.setItem(i, 4, QtWidgets.QTableWidgetItem(str(b[i])))
#             self.tableWidget.setItem(i, 5, QtWidgets.QTableWidgetItem(str(c[i])))
#             self.tableWidget.setItem(i, 6, QtWidgets.QTableWidgetItem(str(d[i])))
#             i+=1
#         i=0
#         while i<=N:
#             k=1
#             if(xAdd[i] > x[k]):
#                 k+=1
#             self.tableWidget_2.setItem(i, 0, QtWidgets.QTableWidgetItem(str(i)))
#             self.tableWidget_2.setItem(i, 1, QtWidgets.QTableWidgetItem(str(xAdd[i])))
#             self.tableWidget_2.setItem(i, 2, QtWidgets.QTableWidgetItem(str(f(xAdd[i]))))
#             self.tableWidget_2.setItem(i, 3, QtWidgets.QTableWidgetItem(str(S(x[k], xAdd[i],a[k],b[k],c[k],d[k]))))
#             self.tableWidget_2.setItem(i, 4, QtWidgets.QTableWidgetItem(str(f(xAdd[i]) - S(x[k], xAdd[i],a[k],b[k],c[k],d[k]))))
#             self.tableWidget_2.setItem(i, 5, QtWidgets.QTableWidgetItem(str(df(xAdd[i]))))
#             self.tableWidget_2.setItem(i, 6, QtWidgets.QTableWidgetItem(str(df(xAdd[i]) - dS(x[k], xAdd[i],a[k],b[k],c[k],d[k]))))

#             self.tableWidget_3.setItem(i, 0, QtWidgets.QTableWidgetItem(str(i)))
#             self.tableWidget_3.setItem(i, 1, QtWidgets.QTableWidgetItem(str(xAdd[i])))
#             self.tableWidget_3.setItem(i, 2, QtWidgets.QTableWidgetItem(str(ddf(xAdd[i]))))
#             self.tableWidget_3.setItem(i, 3, QtWidgets.QTableWidgetItem(str(ddS(x[k], xAdd[i],a[k],b[k],c[k],d[k]))))
#             self.tableWidget_3.setItem(i, 4, QtWidgets.QTableWidgetItem(str(ddf(xAdd[i])-ddS(x[k], xAdd[i],a[k],b[k],c[k],d[k]))))

#             i+=1
#         self.tableWidget.resizeColumnsToContents()
#         self.tableWidget_2.resizeColumnsToContents()
#         self.tableWidget_3.resizeColumnsToContents()
#         secwin.label_2.setText(str(n))
#         secwin.label_2.setText(str(N))
        
#         #spravka######################################################
#         t=0
#         maxFS=0
#         maxdFdS=0
#         maxddFddS=0
#         fs_x=0
#         dfds_x=0
#         ddfdds_x=0
#         while t<=N:
#             k=1
#             if(xAdd[t] > x[k]):
#                 k+=1
#             if( maxFS>abs(f(xAdd[t])-S(x[k], xAdd[t], a[k],b[k],c[k],d[k])) ):
#                 maxFS=abs(f(xAdd[t])-S(x[k], xAdd[t], a[k],b[k],c[k],d[k]))
#                 fs_x=xAdd[t]

#             if( maxdFdS>abs(df(xAdd[t])-dS(x[k], xAdd[t], a[k],b[k],c[k],d[k])) ):
#                 maxdFdS=abs(df(xAdd[t])-dS(x[k], xAdd[t], a[k],b[k],c[k],d[k]))
#                 dfds_x=xAdd[t]

#             if( maxddFddS>abs(ddf(xAdd[t])-ddS(x[k], xAdd[t], a[k],b[k],c[k],d[k])) ):
#                 maxddFddS=abs(ddf(xAdd[t])-ddS(x[k], xAdd[t], a[k],b[k],c[k],d[k]))
#                 ddfdds_x=xAdd[t]
#             t+=1
        
#         secwin.textBrowser.setText("max|F(xj)− S(xj)| = " + str(maxFS) + " при x = " + str(fs_x))
#         secwin.textBrowser_2.setText("max|F'(xj)− S'(xj)| = " + str(maxdFdS)+ " при x = " + str(dfds_x))
#         secwin.textBrowser_3.setText("max|F''(xj)− S''(xj)| = " + str(maxddFddS)+ " при x = " + str(ddfdds_x))
#         ########################################################################################
#         #ГРАФИКИ####################################
#         # ax = self.figure.add_subplot(221)
#         # axd = self.figure.add_subplot(223)
#         # axdd = self.figure.add_subplot(122)

#         ax = self.figure.add_subplot(111)
#         ax.axis([-1.5, 4.5, -1, 3])
#         # axd.axis([-2, 5, -1, 3])
#         # axdd.axis([-2, 5, -1, 3])


#         Spl=np.zeros(N+1)
#         t=0
#         while t<=N:
#             k=1
#             if(xAdd[t] > x[k]):
#                 k+=1
#             Spl[t]=a[k]+b[k]*(xAdd[t] - x[k])+c[k]/2.0*(xAdd[t] - x[k])**2+d[k]/6.0 * (xAdd[t] - x[k])**3
#             t+=1
#         #основной график
#         plt.subplot(111)
#         ax.plot(x, y, '-r')
#         ax.plot(xAdd, Spl, '-b')
        
#         #просто сплайн
        
#         #plt.plot(Spl)

#         #plt.show()

#         #x, xs, y, S, dS, ddS, f,df,ddf = [], [],[],[],[],[],[],[],[]
#         # i=0
#         # while i <= N:
#         #     k=1
#         #     if(xAdd[i] > xList[k]):
#         #         k+=1            
#         #     xs.append(xAdd[i])
#         #     S.append(S(xList[k], xAdd[i], a[k],b[k],c[k],d[k]))
#         #     dS.append(dS(xList[k], xAdd[i], a[k],b[k],c[k],d[k]))
#         #     ddS.append(ddS(xList[k], xAdd[i], a[k],b[k],c[k],d[k]))
#         #     i += 1
#         # i=0
#         # while i <= n:
#         #     x.append(xList[i])
#         #     f.append(yList[i])
#         #     df.append(df(xList[i]))
#         #     ddf.append(ddf(xList[i]))
#         #     i += 1

#         #plt.subplot(111)
        




#         # ax.plot(xAdd, S, '-r')
#         # ax.plot(x, f, '-b')
#         #plt.plot(xdop, S, '-r')

#         # axd.plot(Sxlist, dSlist, '-r')
#         # axd.plot(xlist, dflist, '-b')

#         # axdd.plot(Sxlist, ddSlist, '-r')
#         # axdd.plot(xlist, ddflist, '-b')
#         #ax_2.plot(xlist, u2list, '-b')

#         # ax.legend()
#         # axd.legend()
#         # axdd.legend()
#         ax.grid(True)
#         # axd.grid(True)
#         # axdd.grid(True)
#         self.canvas.draw()










