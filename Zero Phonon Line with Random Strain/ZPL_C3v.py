import numpy as np
import ctypes
import matplotlib.pyplot as plt
#import sys
from scipy.integrate import nquad


def nm2mev(data, dlambda):
   spectrum = data
   spectrum[:,0] += dlambda
   wavelength = spectrum[:,0]
   lambda0 = spectrum[0,0]
   w0 = 1e7/lambda0
   i=0
   for w in wavelength:
     spectrum[i,0] = w0 + (1/w - 1/lambda0)*1e7; #cm-1
     i += 1

   spectrum[:,0] = spectrum[:,0] * 1.23981e-4 * 1e3; #meV
   return spectrum



def gen_freq(xc):
   freq1 = np.arange(-1,1.1,0.1)
   freq = []
   freq = [-5, -4.5, -4, -3.5, -3, -2.5, -2.25, -2, -1.75, -1.5, -1.25]
   freq = np.append(freq,freq1)
   freq = np.append(freq,[1.25, 1.5, 1.75, 2, 2.25, 2.5, 3, 3.5, 4, 4.5, 5])
   freq = freq  + xc
   return freq

   
   
def integrate(freq,delta,gamma,v1,v2,xc):
   y = []
   for x in freq:
      print(x)
      int = nquad(func, [[0, 1],[0, 1],[0, 1],[0, 1]], args=(x,delta,gamma,v1,v2,xc))
      y = np.append(y,int[0])
   print(freq.shape,y.shape)
   return y


if __name__ == "__main__":

   print('loading data')
   filename = 'E6P1_532'
   data = np.genfromtxt(filename+'.csv', delimiter=',')
   dlambda = 1041.6893 
   #dlambda is needed if the spectrum has been translated during preprocessing
   #otherwise, put dlambda = 0
   spectrum = nm2mev(data, dlambda)
   
   #Simulation parameters
   xc = 1190.1917
   freq = gen_freq(xc)
   #the numbers are from Rogers Doherty 2015 New J Phys:
   B = -1.23
   C = -0.69
   c11 = 1076
   c12 = 125
   c44 = 576
   #the following is consistent with the notations in Davies 1979 J Phys C :
   v1 = (c11-c12)*B + c44*C #in meV
   v2 = np.sqrt(2)*(c11-c12)*B - c44*C/np.sqrt(2) #in meV
   #the fitted width parameters for E6P1_532 :
   delta = 2.9e-1
   gamma = 1.7e-4
   
   
   print('loading first dll')
   lib = ctypes.CDLL('lineshape.so')
   func = lib.h
   func.restype = ctypes.c_double
   func.argtypes = (ctypes.c_int, ctypes.c_double)

   print('calculating lineshape')
   y = integrate(freq,delta,gamma,v1,v2,xc)
   
   output = np.column_stack((freq, y))
   np.savetxt('simul_'+filename+'.csv',output,delimiter=",")

   plt.plot(freq, y/max(y), color='r', label='sim')
   plt.hold('on')
   plt.plot(spectrum[:,0], spectrum[:,1]/max(spectrum[:,1]), color='b', label='exp')
   plt.xlabel('Energy, meV')
   plt.ylabel('Normalized intensity')
   plt.show()
