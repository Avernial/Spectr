import numpy as np
import pyqtgraph as pg
import scipy.signal as sig
from pyqtgraph import exporters
import matplotlib.pyplot as plt
#
##############################
#
name_filtr = 'Filtr №1, Scan speed: 128 nm/min'
start_nm = 400
speed_nm = 128
s = 'Scan speed: 128'
name_file = 'test1'

##############################
with open('filtr1-raw.dat', 'rb') as f:
    dtype = np.dtype([('signals','uint8'),('value', 'uint16'),('time','uint32')])
    data = np.fromfile(f,dtype) 
#################
ch1 = (data['signals'] & 0b00000100) >> 2
ch2 = (data['signals'] & 0b00000010) >> 1
time = data['time']
val = data['value']

def beginning_ch(signal):
    m = []
    old_item = False
    for indx, item in enumerate(signal):          
        if item == 1 and old_item == 0:
             m.append(indx)
             break
    return m   
bch = beginning_ch(ch1)

ch1 = list(ch1[bch[0]:])
ch2 = list(ch2[bch[0]:])
val = list(val[bch[0]:])
time = list(time[bch[0]:])

def range_ch(signal):
    m = []
    n = []
    old_item = False
    old_indx = False
    i = 0
    while i <= len(signal):
        if signal[i] == 1 and old_item == 0:
            m.append(old_indx)
            i+=260
            j = i
            while True:
                j+=-1
                if signal[j] == 1:
                    n.append(j)
                    break
        old_item = signal[i] 
        old_indx = i
        i+=1    
        if i+260 > len(signal):
            break
        
    return m, n

rch1 = range_ch(ch1)
rch2 = range_ch(ch2)

def division1(x, y, m, n):
    v = []
    t = []       
    for sl in zip(m, n):
            v.append(x[sl[0]:sl[1]])
            t.append(y[sl[0]:sl[1]])
    return v, t

div1 = division1(val, time, rch1[0], rch1[1]) 
div2 = division1(val, time, rch2[0], rch2[1])

ch1 = div1[0]
ch2 = div2[0]

time1 = div1[1]
time2 = div2[1]

with open('2obr-raw-ch1.txt', 'wb') as f:
    np.save(f, np.asarray(ch1))

with open('2obr-raw-ch2.txt', 'wb') as f:    
    np.save(f, np.asarray(ch2))
    
with open('2obr-raw-time1.txt', 'wb') as f:    
    np.save(f, np.asarray(time1))
    
with open('2obr-raw-time2.txt', 'wb') as f:    
    np.save(f, np.asarray(time2)) 

def time_to_nm(time):
    m = []
    for i in time:
        mean = np.mean(i)
        nm = (mean/1000000)*(speed_nm/60)
        m.append(nm)
    return m

nm1 = time_to_nm(div1[1])
nm2 = time_to_nm(div2[1])

nm1 = np.asarray(nm1)
nm2 = np.asarray(nm2)
nm1 = (nm1/2)+start_nm
nm2 = (nm2/2)+start_nm
nm = (nm1+nm2)/2
nm = nm[1:]
with open('2obr-raw-nm.txt', 'wb') as f:    
    np.save(f, np.asarray(nm)) 

data1 = np.load("2obr-raw-ch1.txt")
data2 = np.load("2obr-raw-ch2.txt")
time1 = np.load("2obr-raw-time1.txt")
time2 = np.load("2obr-raw-time2.txt")

def savgol_filtr(data):
    sr = []  
    for i in range(len(data)):
        arr = data[i]
        arr = arr[:int(len(arr)*0.7)]
        sarr = np.sort(arr)
        sarr1 = sarr[:int(len(arr)*0.9)]
        filt_sarr = sig.savgol_filter(sarr1, 15, 1, deriv = 0, delta = 1.0, axis = -1, mode = 'interp', cval = 0.0)
        filt_diff = sig.savgol_filter(np.diff(filt_sarr)*10, 15, 1, deriv = 0, delta = 1.0, axis = -1, mode='interp', cval=0.0)
        ind = np.where(filt_diff == np.max(filt_diff))
        plato = sarr[ind[0][0]:int(len(arr)*0.92)]
#        plato = plato[30:]        
        sr.append(np.mean(plato))
    
    return sr
        
U1 = savgol_filtr(data1)
U2 = savgol_filtr(data2)

u1 = np.asarray(U1)
u2 = np.asarray(U2)
T = u2/u1
T = T[1:]*100
filt_T = sig.savgol_filter(T, 35, 1, deriv = 0, delta = 1.0, axis = -1, mode = 'interp', cval = 0.0)
filt_T2 = sig.savgol_filter(T, 15, 1, deriv = 0, delta = 1.0, axis = -1, mode = 'interp', cval = 0.0)

with open('2obr-raw_savgalov_filt.txt', 'wb') as f:
    np.save(f, np.asarray(T))
    
nm = np.load("2obr-raw-nm.txt")    
T = np.load("2obr-raw_savgalov_filt.txt")    

np.savez('Filtr2.npz', Wavelength = nm, T = T)
data = np.load('Filtr2.npz')
#filt_T = sig.savgol_filter(data['Transmittance'], 21, 1, deriv = 0, delta = 1.0, axis = -1, mode = 'interp', cval = 0.0)
#filt_T2 = sig.savgol_filter(data['Transmittance'], 21, 1, deriv = 0, delta = 1.0, axis = -1, mode = 'interp', cval = 0.0)
#plt.plot(data['Nanometers'], data['Transmittance'])
#plt.savefig('C:\\Users\\FRostNOva\\Desktop\\{}'.format(name_filtr))

win = pg.GraphicsWindow(title="Spectr")

win.resize(1000, 600)
p1 = win.addPlot(name = "Plot1", title = name_filtr)

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
#p1.plot(data['Nanometers'],filt_T, pen = 'r')
#p1.plot(data['Nanometers'], filt_T2*100, pen = 'k')
p1.plot(data['Wavelength'], data['T'], pen = 'k')
p1.setLabel('left', "T, %")
p1.setLabel('bottom', "Wavelength, nm")
p1.showGrid(x = True, y = True, alpha = 0.7)
p1.setYRange(0, 100)
pg.QtGui.QApplication.exec_()


#exporter = pg.exporters.ImageExporter(p1)
#exporter.parameters()['width'] = 1920
#exporter.parameters()['height'] = 1080
#params = exporter.parameters()
#params.param("width").setValue(1920, blockSignal = exporter.widthChanged)
#params.param("height").setValue(1080, blockSignal = exporter.heightChanged)
#exporter.export('fileName1.png')