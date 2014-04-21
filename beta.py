from numpy import genfromtxt, array,linspace, exp, log,sqrt
from uncertainties import  ufloat
from matplotlib.pyplot import *

from Tools import lin_reg , make_LaTeX_table
data = genfromtxt("beta_data1.csv", delimiter=',', unpack="true")

N_3 = data[1]
d_3 = data[0]*1e-6
N_0 = 1556
t_0= 900.0
t_3 = 60.0

def fehler(N,N_0,t,t_0):
    return sqrt(N/t**2 + N_0/t_0**2)


A_0= N_0/t_0
A_3 = N_3/t_3
A_3_err = array([fehler(i,N_0,t_3,t_0) for i in N_3])




m_1,b_1 = lin_reg(d_3[:12],log(A_3[:12]))
m_2,b_2 = lin_reg(d_3[13:], log(A_3[13:]))

x_1 = linspace(0,1e-3)
x_2 = linspace(0,2500e-6)



xlabel("Schichtdicke [$\mu m$]")
ylabel("$A$ in [$s^{-1}$]")

errorbar(d_3*1e6, A_3, A_3_err,0,'.')
plot(x_1*1e6,exp( m_1.n*x_1+b_1.n))
plot(x_2*1e6, exp(m_2.n*x_2+b_2.n))
yscale("log")
show()
