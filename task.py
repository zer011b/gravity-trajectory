#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys
from PyQt5.QtWidgets import QWidget, QApplication, QPushButton, QLineEdit, QCheckBox
from PyQt5.QtGui import QPainter, QColor, QPen
from PyQt5.QtCore import Qt, QPointF
from math import exp

# const
g=9.8

# lists with values
x=[]
y=[]
vx=[]
vy=[]

x_exact=[]
y_exact=[]
vx_exact=[]
vy_exact=[]

# default values for global variables
N=50
dt=0.1
k=1.0
m=1.0
x0=0.0
y0=10.0
vx0=1.0
vy0=0.0

# parameters
useK=1

# window size
window0x=300
window0y=300
window_sizex=900
window_sizey=500

# some window modifiers
diffx=50
diffy=50

def calculate_k1(N, dt, k, m, x0, y0, vx0, vy0):
  global x, y, vx, vy

  x = [0] * N
  y = [0] * N
  vx = [0] * N
  vy = [0] * N

  x[0] = x0
  y[0] = y0
  vx[0] = vx0
  vy[0] = vy0

  for i in range(1,N):
    vx[i] = vx[i-1] - k/m * vx[i-1] * dt
    vy[i] = vy[i-1] - (g + k/m * vy[i-1]) * dt

  for i in range(1,N):
    x[i] = x[i-1] + vx[i] * dt
    y[i] = y[i-1] + vy[i] * dt


def calculate_k1_exact(N, dt, k, m, x0, y0, vx0, vy0):
  global x_exact, y_exact, vx_exact, vy_exact

  x_exact = [0] * N
  y_exact = [0] * N
  vx_exact = [0] * N
  vy_exact = [0] * N

  for i in range(0,N):
    vx_exact[i] = vx0 * exp (-k*i*dt/m)
    vy_exact[i] = vy0 * exp (-k*i*dt/m) - g*m/k* (1 - exp (-k*i*dt/m))
    x_exact[i] = x0 + vx0 * m / k * (1 - exp (-k*i*dt/m))
    y_exact[i] = y0 + m / k * ((vy0 + m*g/k) * (1 - exp (-k*i*dt/m)) - g*i*dt)


def calculate_k0(N, dt, x0, y0, vx0, vy0):
  global x, y, vx, vy

  x = [0] * N
  y = [0] * N
  vx = [0] * N
  vy = [0] * N

  x[0] = x0
  y[0] = y0
  vx[0] = vx0
  vy[0] = vy0

  for i in range(1,N):
    vx[i] = vx[i-1]
    vy[i] = vy[i-1] - g * dt

  for i in range(1,N):
    x[i] = x[i-1] + vx[i] * dt
    y[i] = y[i-1] + vy[i] * dt

def calculate_k0_exact(N, dt, x0, y0, vx0, vy0):
  global x_exact, y_exact, vx_exact, vy_exact

  x_exact = [0] * N
  y_exact = [0] * N
  vx_exact = [0] * N
  vy_exact = [0] * N

  for i in range(0,N):
    vx_exact[i] = vx0
    vy_exact[i] = vy0 - g * i *dt
    x_exact[i] = x0 + vx0 * i * dt
    y_exact[i] = y0 + vy0 * i * dt - g * i * i * dt * dt /2


def calculate(N, dt, k, m, x0, y0, vx0, vy0):
  if useK == 1:
    calculate_k1(N, dt, k, m, x0, y0, vx0, vy0)
  else:
    calculate_k0(N, dt, x0, y0, vx0, vy0)

def calculate_exact(N, dt, k, m, x0, y0, vx0, vy0):
  if useK == 1:
    calculate_k1_exact(N, dt, k, m, x0, y0, vx0, vy0)
  else:
    calculate_k0_exact(N, dt, x0, y0, vx0, vy0)



class Example(QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()

    def button_click (self):
      global N, dt, k, m, x0, y0, vx0, vy0, useK
      N=int(self.le1.text())
      dt=float(self.le2.text())
      k=float(self.le3.text())
      m=float(self.le4.text())
      x0=float(self.le5.text())
      y0=float(self.le6.text())
      vx0=float(self.le7.text())
      vy0=float(self.le8.text())

      if k == 0.0:
        print(useK)
        useK = 0
        self.cb.setChecked(Qt.Unchecked)
        self.update()

      calculate (N, dt, k, m, x0, y0, vx0, vy0)
      calculate_exact (N, dt, k, m, x0, y0, vx0, vy0)
      self.update()

    def changeUseK (self, state):
      global useK
      
      if state == Qt.Checked:
        useK=1
      else:
        useK=0

    def initUI(self):

        self.setGeometry(window0x, window0y, window_sizex, window_sizey)
        self.setWindowTitle('Gravity')

        self.cb = QCheckBox('Use K', self)
        self.cb.move(330, 30)
        self.cb.toggle ()
        self.cb.stateChanged.connect(self.changeUseK)

        self.pb = QPushButton("calculate", self)
        self.pb.move (20, 10)

        self.le1 = QLineEdit(str(N), self)
        self.le1.move (140, 10)
        self.le1.setFixedWidth (50)

        self.le2 = QLineEdit(str(dt), self)
        self.le2.move (230, 10)
        self.le2.setFixedWidth (50)

        self.le3 = QLineEdit(str(k), self)
        self.le3.move (330, 10)
        self.le3.setFixedWidth (50)

        self.le4 = QLineEdit(str(m), self)
        self.le4.move (430, 10)
        self.le4.setFixedWidth (50)

        self.le5 = QLineEdit(str(x0), self)
        self.le5.move (530, 10)
        self.le5.setFixedWidth (50)

        self.le6 = QLineEdit(str(y0), self)
        self.le6.move (630, 10)
        self.le6.setFixedWidth (50)

        self.le7 = QLineEdit(str(vx0), self)
        self.le7.move (730, 10)
        self.le7.setFixedWidth (50)

        self.le8 = QLineEdit(str(vy0), self)
        self.le8.move (830, 10)
        self.le8.setFixedWidth (50)

        self.pb.clicked.connect(self.button_click)

        self.show()


    def paintEvent(self, e):

        qp = QPainter()
        qp.begin(self)
        self.drawLines(qp)
        qp.end()

    def drawLines(self, qp):

        maxx=x[N-1]
        if maxx == 0:
          maxx = 1

        shiftx=((window_sizex-2*diffx)/maxx)
        # todo: make 4*max(y[0], y_exact[0]) or even max(y, y_exact)
        shifty=((window_sizey-2*diffy)/(4* y[0]))

	# draw some string
        pen = QPen(Qt.black, 1, Qt.SolidLine)
        qp.setPen(pen)

        qp.drawText (QPointF(diffx+0+5, window_sizey/2-5), "O")
        qp.drawText (QPointF(diffx-10,diffy), "Y")

        # make
        qp.drawText (QPointF(window_sizex-diffx, window_sizey/2-5), "X")
        qp.drawText (QPointF(120, 28), "N=")
        qp.drawText (QPointF(206, 28), "dt=")
        qp.drawText (QPointF(310, 28), "k=")
        qp.drawText (QPointF(410, 28), "m=")
        qp.drawText (QPointF(510, 28), "x0=")
        qp.drawText (QPointF(610, 28), "y0=")
        qp.drawText (QPointF(710, 28), "vx0=")
        qp.drawText (QPointF(810, 28), "vy0=")

        # todo: fix this
        qp.drawText (QPointF(diffx+5,diffy), str(2*y[0]))
        qp.drawText (QPointF(diffx+5,window_sizey - diffy), str(-2*y[0]))
        qp.drawText (QPointF(diffx+shiftx*x[N-1] - 5, window_sizey/2+20), str(x[N-1]))

        pen = QPen(Qt.black, 2, Qt.DashLine)
        qp.setPen(pen)
        qp.drawLine(diffx+0, window_sizey/2, diffx+shiftx*x[N-1], window_sizey/2)
        qp.drawLine(diffx+0, window_sizey-diffy-0, diffx+0, diffy+0)

        pen = QPen(Qt.red, 1, Qt.SolidLine)
        qp.setPen(pen)

        for i in range(1,N):
          qp.drawLine(diffx+shiftx*x[i-1], window_sizey/2-shifty*y[i-1], diffx+shiftx*x[i], window_sizey/2-shifty*y[i])

        pen = QPen(Qt.green, 1, Qt.SolidLine)
        qp.setPen(pen)

        for i in range(1,N):
          qp.drawLine(diffx+shiftx*x_exact[i-1], window_sizey/2-shifty*y_exact[i-1], diffx+shiftx*x_exact[i], window_sizey/2-shifty*y_exact[i])


if __name__ == '__main__':

    app = QApplication(sys.argv)
    calculate(N, dt, k, m, x0, y0, vx0, vy0)
    calculate_exact(N, dt, k, m, x0, y0, vx0, vy0)
    ex = Example()
    sys.exit(app.exec_())
