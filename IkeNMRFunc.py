import numpy as np
import scipy as sp
import scipy.linalg
import matplotlib.pyplot as plt
from math import *
import cmath
def fact(n):
    if n == 0:
        return 1
    else:
        return n * fact(n-1)

def sqrt(n): return n**(1/2)
def exp(n): return cmath.exp(n)
def expm(n): return sp.linalg.expm(n)
def trace(n): return np.trace(n)

def Wig(theta):
    mp_a = np.linspace(-I,I,int(2*I+1))
    m_a  = np.linspace(-I,I,int(2*I+1))
    s_a = np.linspace(0,2*I,int(2*I+1))
    d = np.zeros((int(2*I+1),int(2*I+1)))+1j*np.zeros((int(2*I+1),int(2*I+1)))
    for mp in mp_a:
        for m in m_a:
            sum = 0
            for s in s_a:
                if I+m-s >=0 and s >=0 and mp-m+s >=0 and I-mp-s >=0:
                    B = ((-1)**(s)*sqrt(fact(I+mp)*fact(I-mp)*fact(I+m)*fact(I-m)))/(fact(I+m-s)*fact(s)*fact(mp-m+s)*fact(I-mp-s))
                    C = (cos(theta/2))**(2*I+m-mp-(2*s))
                    D = (sin(theta/2))**(mp-m+(2*s))
                    sum = sum+B*exp(-1*1j*np.pi/2*mp)*C*D*exp(1j*np.pi/2*mp)
            d[int(mp+I),int(m+I)] = sum
    return d

def LVN(R,M): return np.matmul(np.matmul(R,M),np.linalg.inv(R))
def SAND(R1,M,R2): return np.matmul(np.matmul(R,M),np.linalg.inv(R))
#Generating Pauli Matrices for a given spin state I
def PauliMatrices(I):
    j_a = np.linspace(-I,I-1,int(2*I))
    k_a = np.linspace(-I,I,int(2*I+1))
    N = int(2*I+1)
    Splus  = np.zeros((N,N))
    Sminus = np.zeros((N,N))
    Sz = np.zeros((N,N))
    for j in j_a:
        m_s = j
        Splus[int(j+I),int(j+I+1)] = (I*(I+1)-m_s*(m_s+1))**(1/2)
        Sminus[int(j+I+1),int(j+I)] = (I*(I+1)-(m_s+1)*m_s)**(1/2)

    Sx = (Splus+Sminus)/2
    Sy = -1j*((Splus-Sminus))/2

    for k in k_a:
        Sz[(int(-k+I),int(-k+I))] = k
    return [Sx,Sy,Sz,Splus,Sminus]
