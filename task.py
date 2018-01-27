#!/usr/bin/python3
#-*- coding: utf-8 -*-

# import system libs
# подключение системных библиотек
import sys
from PyQt5.QtWidgets import QWidget, QApplication, QPushButton, QLineEdit, QCheckBox, QComboBox
from PyQt5.QtGui import QPainter, QColor, QPen
from PyQt5.QtCore import Qt, QPointF
from math import exp, sqrt

# constants
# константы
g = 9.8

# default values for global variables
# глобальные переменные с их значениями по умолчанию

# number of points
# число точек
N = 50
# numerical step
# разностный шаг по времени
dt = 0.1
# coefficient of air resistance
# коэффициент сопротивления воздуха
k = 1.0
# mass
# масса тела
m = 1.0
# default coordinate
# начальные координаты
x0 = 0.0
y0 = 10.0
# default speed
# начальная скорость
vx0 = 1.0
vy0 = 0.0

# order of accuracy of numerical method
# порядок точности численного метода
orderOfAccuracy = 1

# flag, whether to use air resistance
# флаг, использовать ли сопротивление воздуха
useAirResistance = True

# default window position
# начальная позиция окна
window0x = 50
window0y = 50

# default window size
# начальные размеры окна
window_sizex = 950
window_sizey = 700

# lists with numerical values
# списки со значениями координат, полученных численно
x = []
y = []
# списки со значениями скоростей, полученных численно
vx = []
vy = []

# lists with exact values
# списки со значениями координат, полученных из точного решения
x_exact = []
y_exact = []
# списки со значениями скоростей, полученных из точного решения
vx_exact = []
vy_exact = []

# Норма разности точного и численного решений
norm = 0

# some window modifiers
diffx=50
diffy=100

# функция для расчета нормы разности точного и численного решения
def calculate_norm ():
  global x, y, x_exact, y_exact, norm

  norm = 0

  for i in range (1,N):
    norm += (y[i] - y_exact[i]) ** 2

  norm = sqrt(norm)


# function to calculate numerically with air resistance
# функция для численного расчета с учетом сопротивления воздуха
def calculate_k1 (N, dt, k, m, x0, y0, vx0, vy0, order):
  global x, y, vx, vy

  x = [0] * N
  y = [0] * N
  vx = [0] * N
  vy = [0] * N

  x[0] = x0
  y[0] = y0
  vx[0] = vx0
  vy[0] = vy0

  for i in range (1,N):
    if order == 1:
      # метод Эйлера первого порялка точности
      vx[i] = vx[i-1] - k/m * vx[i-1] * dt
      vy[i] = vy[i-1] - (g + k/m * vy[i-1]) * dt
    elif order == 2:
      # метод Рунге-Кутта второго порядка точности
      vx[i] = vx[i-1] * (1 - k/m * dt + 0.5 * (k/m * dt) ** 2)
      vy[i] = (k/m * g * (dt)**2 / 2 - dt * g) + vy[i-1] * (1 - k/m * dt + 0.5 * (k/m * dt) ** 2)

  for i in range (1,N):
    if order == 1:
      # метод Эйлера первого порялка точности
      x[i] = x[i-1] + vx[i] * dt
      y[i] = y[i-1] + vy[i] * dt
    elif order == 2:
      # метод Рунге-Кутта второго порядка точности
      x[i] = x[i-1] + dt/2 * (vx[i] + vx[i-1])
      y[i] = y[i-1] + dt/2 * (vy[i] + vy[i-1])

# function to calculate exact with air resistance
# функция для точного расчета с учетом сопротивления воздуха
def calculate_k1_exact (N, dt, k, m, x0, y0, vx0, vy0):
  global x_exact, y_exact, vx_exact, vy_exact

  x_exact = [0] * N
  y_exact = [0] * N
  vx_exact = [0] * N
  vy_exact = [0] * N

  for i in range (0,N):
    vx_exact[i] = vx0 * exp (-k*i*dt/m)
    vy_exact[i] = vy0 * exp (-k*i*dt/m) - g*m/k* (1 - exp (-k*i*dt/m))
    x_exact[i] = x0 + vx0 * m/k * (1 - exp (-k*i*dt/m))
    y_exact[i] = y0 + m / k * ((vy0 + m*g/k) * (1 - exp (-k*i*dt/m)) - g*i*dt)

# function to calculate numerically without air resistance
# функция для численного расчета без учета сопротивления воздуха
def calculate_k0 (N, dt, x0, y0, vx0, vy0, order):
  global x, y, vx, vy

  x = [0] * N
  y = [0] * N
  vx = [0] * N
  vy = [0] * N

  x[0] = x0
  y[0] = y0
  vx[0] = vx0
  vy[0] = vy0

  for i in range (1,N):
    vx[i] = vx[i-1]
    vy[i] = vy[i-1] - g * dt

  for i in range (1,N):
    if order == 1:
      # метод Эйлера первого порялка точности
      x[i] = x[i-1] + vx[i] * dt
      y[i] = y[i-1] + vy[i] * dt
    elif order == 2:
      # метод Рунге-Кутта второго порядка точности
      x[i] = x[i-1] + dt/2 * (vx[i] + vx[i-1])
      y[i] = y[i-1] + dt/2 * (vy[i] + vy[i-1])

# function to calculate exact without air resistance
# функция для точного расчета без учета сопротивления воздуха
def calculate_k0_exact (N, dt, x0, y0, vx0, vy0):
  global x_exact, y_exact, vx_exact, vy_exact

  x_exact = [0] * N
  y_exact = [0] * N
  vx_exact = [0] * N
  vy_exact = [0] * N

  for i in range (0,N):
    vx_exact[i] = vx0
    vy_exact[i] = vy0 - g * i *dt
    x_exact[i] = x0 + vx0 * i * dt
    y_exact[i] = y0 + vy0 * i * dt - g * i * i * dt * dt /2

# function to calculate numerically
# функция для численного расчета
def calculate (N, dt, k, m, x0, y0, vx0, vy0, order):
  if useAirResistance:
    calculate_k1 (N, dt, k, m, x0, y0, vx0, vy0, order)
  else:
    calculate_k0 (N, dt, x0, y0, vx0, vy0, order)

# function to calculate exact
# функция для точного расчета
def calculate_exact (N, dt, k, m, x0, y0, vx0, vy0):
  if useAirResistance:
    calculate_k1_exact (N, dt, k, m, x0, y0, vx0, vy0)
  else:
    calculate_k0_exact (N, dt, x0, y0, vx0, vy0)


# QT widget to draw GUI
# Виджет для рисования графического интерфейса
class TaskWidget (QWidget):

    # constructor of TaskWidget
    # конструктор объектов класса TaskWidget
    def __init__(self):
        # call constructor on parent object
        # функция super() возвращает родительский объект, и мы вызываем его конструктор
        super().__init__()

        # create GUI
        # вызываем функцию, создающую графический интерфейс
        self.initUI()

    # функция для обработки выбора численного метода
    def change_order (self, index):
      global orderOfAccuracy
      orderOfAccuracy = index + 1


    # button click handler
    # функция для обработки нажатия кнопки
    def button_click (self):
      global N, dt, k, m, x0, y0, vx0, vy0, useAirResistance

      # get values from editor windows
      # получаем введенные значения
      N=int(self.le1.text())
      dt=float(self.le2.text())
      k=float(self.le3.text())
      m=float(self.le4.text())
      x0=float(self.le5.text())
      y0=float(self.le6.text())
      vx0=float(self.le7.text())
      vy0=float(self.le8.text())

      # if coefficient k is set to zero, use non-resistance mode
      # переключаем режим, если коэффициент k был задан равным 0
      if k == 0.0:
        useAirResistance = False
        self.cb.setChecked(Qt.Unchecked)
        self.update()

      # calculate numerical and exact solutions
      # вычисляем численное и точное решения
      calculate (N, dt, k, m, x0, y0, vx0, vy0, orderOfAccuracy)
      calculate_exact (N, dt, k, m, x0, y0, vx0, vy0)
      calculate_norm()

      # update GUI
      # обновляем интерфейс
      self.update()


    # toggle handler
    # функция для обработки переключения режима расчета с/без сопротивлением воздуха
    def changeUseAirResistance (self, state):
      global useAirResistance

      # set mode according to toggle
      # задаем режим согласно переключателю в графическом интерфейсе
      if state == Qt.Checked:
        useAirResistance=True
      else:
        useAirResistance=False


    # create GUI
    # функция для создания графического интерфейса
    def initUI(self):
        # geometry of window
        # задание геометрии окна и заголовка
        self.setGeometry(window0x, window0y, window_sizex, window_sizey)
        #self.setWindowTitle('Gravity')
        self.setWindowTitle('Движение тела в поле силы тяжести')

        # toggle for resistance
        # задание переключателя режима с/без сопротивления
        #self.cb = QCheckBox('Use air resistance', self)
        self.cb = QCheckBox('Использовать сопротивление воздуха', self)
        self.cb.move(330, 40)
        self.cb.toggle ()
        self.cb.stateChanged.connect(self.changeUseAirResistance)

        # button for calculation
        # кнопка для начала расчета
        #self.pb = QPushButton("calculate", self)
        self.pb = QPushButton("Расчет", self)
        self.pb.move (20, 10)
        self.pb.clicked.connect(self.button_click)

        self.сb = QComboBox(self)
        self.сb.addItems(['Первый порядок точности', 'Второй порядок точности'])
        self.сb.move (20, 40)
        self.сb.currentIndexChanged.connect(self.change_order)

        # edit fields
        # поля для ввода параметров
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

        # show widget
        # показать виджет
        self.show()


    # convert x coordinate to position
    # функция для преобразования координаты x в позицию на экране
    def xToPos(self,val):
        return diffx + self.shiftx * val


    # convert y coordinate to position
    # функция для преобразования координаты y в позицию на экране
    def yToPos(self,val):
        return diffy + self.shifty * max(y_exact) - self.shifty * val


    # init shift coefficients
    # функция для вычисления коэффициентов для рисования
    def initShifts(self):
        maxx = max(x_exact)
        if maxx == 0:
          maxx = 1

        maxy = max(y_exact)
        miny = min(y_exact)
        if maxy == miny:
          maxy = miny + 1

        self.shiftx = (window_sizex - 2*diffx) / maxx
        self.shifty = (window_sizey - 2*diffy) / (maxy - miny)


    # draw strings
    # функция для рисования имен параметров
    def drawStrings(self, qp):
        pen = QPen(Qt.black, 1, Qt.SolidLine)
        qp.setPen(pen)

        # draw strings
        # рисование имен параметров
        qp.drawText (QPointF(120, 28), "N=")
        qp.drawText (QPointF(206, 28), "dt=")
        qp.drawText (QPointF(310, 28), "k=")
        qp.drawText (QPointF(405, 28), "m=")
        qp.drawText (QPointF(503, 28), "x0=")
        qp.drawText (QPointF(603, 28), "y0=")
        qp.drawText (QPointF(697, 28), "vx0=")
        qp.drawText (QPointF(797, 28), "vy0=")

        pen = QPen(Qt.red, 2, Qt.SolidLine)
        qp.setPen(pen)
        qp.drawText (QPointF(700, 60), "Численное решение")

        pen = QPen(Qt.green, 2, Qt.SolidLine)
        qp.setPen(pen)
        qp.drawText (QPointF(700, 80), "Точное решение")

        pen = QPen(Qt.red, 4, Qt.SolidLine)
        qp.setPen(pen)
        qp.drawText (QPointF(20, 650), "Норма разности точного и численного решений: " + str(norm))


    # paint of widget
    # функция рисования для виджета
    def paintEvent(self, e):
        qp = QPainter()
        qp.begin(self)
        self.drawLines(qp)
        qp.end()


    # paint axes and lines
    # функция рисования осей и решений
    def drawLines(self, qp):
        # calculate shift coefficients
        # посчитаем коэффициент для рисования
        self.initShifts()

        # draw strings
        # рисование имен параметров
        self.drawStrings(qp)

        pen = QPen(Qt.black, 1, Qt.SolidLine)
        qp.setPen(pen)

        # draw axes
        # рисование обозначений осей
        qp.drawText (QPointF(self.xToPos(0) + 5, self.yToPos(0) - 5), "O")
        qp.drawText (QPointF(window_sizex - diffx, self.yToPos(0) - 5), "X")
        qp.drawText (QPointF(self.xToPos(0) - 10, diffy), "Y")

        qp.drawText (QPointF(self.xToPos(0) + 5, diffy), str(max(y_exact)))
        qp.drawText (QPointF(self.xToPos(0) + 5, window_sizey - diffy), str(min(y_exact)))
        qp.drawText (QPointF(self.xToPos(x_exact[N-1]) - 5, self.yToPos(0) + 20), str(x_exact[N-1]))

        # draw axes
        # рисование осей
        pen = QPen(Qt.black, 2, Qt.DashLine)
        qp.setPen(pen)

        qp.drawLine(self.xToPos(0), self.yToPos(0), self.xToPos(x_exact[N-1]), self.yToPos(0))
        qp.drawLine(self.xToPos(0), window_sizey - diffy, self.xToPos(0), diffy)

        # draw numerical solution
        # рисование численного решения
        pen = QPen(Qt.red, 1, Qt.SolidLine)
        qp.setPen(pen)

        for i in range(1,N):
          qp.drawLine(self.xToPos(x[i-1]), self.yToPos(y[i-1]), self.xToPos(x[i]), self.yToPos(y[i]))

        # draw exact solution
        # рисование точного решения
        pen = QPen(Qt.green, 1, Qt.SolidLine)
        qp.setPen(pen)

        for i in range(1,N):
          qp.drawLine(self.xToPos(x_exact[i-1]), self.yToPos(y_exact[i-1]), self.xToPos(x_exact[i]), self.yToPos(y_exact[i]))


# starting point
# начало выполнения программы
if __name__ == '__main__':

    # create application
    # создание объекта приложения
    app = QApplication(sys.argv)

    # calculate with default values
    # вычисление с начальными параметрами
    calculate(N, dt, k, m, x0, y0, vx0, vy0, orderOfAccuracy)
    calculate_exact(N, dt, k, m, x0, y0, vx0, vy0)
    calculate_norm()

    # create TaskWidget object and call its constructor
    # создаем объект TaskWidget и вызываем его конструктор
    ex = TaskWidget ()

    # launch app cycle
    # запуск основного цикла событий
    sys.exit(app.exec_())
