"""
ZPL_C3v (by TB): 

The script simulates the zero phonon line for a center
with C3v symmetry in the presence of electron-deformation
interaction with a random strain field produced by point 
defects in the host lattice.
Based on [Malkin et al. PRB 2012].

The lineshape is uploaded as a shared library.
To be generated from lineshape.c file (e.g. with gcc):
  $ gcc -shared -o lineshape.so lineshape.c
   
"""

import numpy as np
import ctypes
import matplotlib.pyplot as plt
#import sys
from scipy.integrate import nquad


def get_spectrum(data, dlambda):
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



def gen_freq():
   freq1 = np.arange(-1,1.1,0.1)
   freq = []
   freq = [-3, -2.5, -2, -1.5]
   freq = np.append(freq,freq1)
   freq = np.append(freq,[1.5, 2, 2.5, 3])  + 1190.1917
   return freq

   
   
def integrate(freq,G,gamma,v1,v2):
   y = []
   for x in freq:
      print(x)
      int = nquad(func, [[0, 1],[0, 1],[0, 1],[0, 1]], args=(x,G,gamma,v1,v2))
      y = np.append(y,int[0])
   print(freq.shape,y.shape)
   return y


if __name__ == "__main__":
   print('loading data')
   filename = 'E6P1_532'
   data = np.genfromtxt(filename+'.csv', delimiter=',')
   dlambda = 1041.6893 
   #dlambda is needed if the spectrum has been translated during processing
   spectrum = get_spectrum(data, dlambda)
   
   freq = gen_freq()
   
   print('loading dll')
   lib = ctypes.CDLL('lineshape.so')
   func = lib.h
   func.restype = ctypes.c_double
   func.argtypes = (ctypes.c_int, ctypes.c_double)


   #the numbers (in cm-1) are from Rogers Doherty 2015 New J Phys:
   B = -1.23
   C = -0.69
   c11 = 1076
   c12 = 125
   c44 = 576
   #the following is from Davies 1979 J Phys C :
   v1 = (c11-c12)*B + c44*C #in meV
   v2 = np.sqrt(2)*(c11-c12)*B - c44*C/np.sqrt(2) #in meV

   #the fitted width parameters for E6P1_532 :
   G = 2.8e-1
   gamma = 2.35e-4
   
   print('calculating lineshape')
   y = integrate(freq,G,gamma,v1,v2)
   
   output = np.column_stack((freq, y))
   np.savetxt('simul'+filename+'.csv',output,delimiter=",")

   plt.plot(freq, y/max(y), color='r', label='sim')
   plt.hold('on')
   plt.plot(spectrum[:,0], spectrum[:,1]/max(spectrum[:,1]), color='b', label='exp')
   plt.xlabel('Wavenumber, nm')
   plt.ylabel('Normalized intensity')
   plt.show()

