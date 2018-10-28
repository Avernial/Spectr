import numpy as np
import pyqtgraph as pg

with open('filtr3-raw.dat', 'rb') as f:
    dtype = np.dtype([('signals','uint8'),('value', 'uint16'),('time','uint32')])
    data = np.fromfile(f,dtype) 
#######
#data = data[:-100000]
b = (data['signals'] & 0b00000100) >> 2
c = (data['signals'] & 0b00000010) >> 1
#d = (data['signals'] & 0b00000001) >> 0
signal = c + b
time = data['time']
val = data['value']
b1 = b * 200
c1 = c * 200
#d1 = d * 100
    
win = pg.GraphicsWindow(title = "Basic plotting")
win.resize(1000,  600)
p1 = win.addPlot(name = "Plot1", title = "Plot1")
#p2 = win.addPlot(name = "Plot2", title = "Plot2: X linked with Plot1", row = 3, col = 0)
#p2.setXLink(p1)
p3 = win.addPlot(name = "Plot3", title = "Plot3: X linked with Plot1", row = 2, col = 0)
p3.setXLink(p1)

p1.plot(time, val, pen = 'c')
#p2.plot(data['time'], d1, pen = 'y')
p3.plot(data['time'], b1, pen = 'r')
p3.plot(data['time'], c1, pen = 'b')
#p1.plot(val, pen = 'c')
#p2.plot(b1, pen = 'y')
#p3.plot(d1, pen = 'r')
#p3.plot(c1, pen = 'b')
pg.QtGui.QApplication.exec_()



arr = np.arange(len(time))
pg.plot(arr, time)
#
#ti = 250
#t =390
#T = ti + t
#
#arr1 = np.zeros(len(c))
#arr2 = np.zeros(len(d))

#for i in range(50, len(arr1), T):
#    arr1[i:i+ti] = 1
#for i in range(350, len(arr2), T):
#    arr2[i:i+ti] = 1
#arr1 = arr1*100
#arr2 = arr2*100