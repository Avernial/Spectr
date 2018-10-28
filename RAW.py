import serial
import time

ser = serial.Serial(port = '\\.\COM4', baudrate = 500000, bytesize = serial.EIGHTBITS, parity = serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE, timeout=1)
time.sleep(4)
with open('filtr3-raw.dat','wb') as f:
    ser.write(b'n')
    print(ser.read(2))
    for _ in range(240000):
        f.write(ser.read(70))
    ser.write(b'n')
    print(ser.read())
    ser.close()       