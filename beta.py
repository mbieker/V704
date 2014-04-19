from numpy import genfromtxt, array,linspace, exp, log
from uncertainties import *
from matplotlib.pyplot import *

from Tools import lin_reg
data = genfromtxt("beta_data1.csv", delimiter=',', unpack="true")

N = data[1]/60
d = data[0]*1e-6
N_0 = 1556/900.

m_1,b_1 = lin_reg(d[:12],log(N[:12]))
m_2,b_2 = lin_reg(d[13:], log(N[13:]))

x_1 = linspace(0,1e-3)
x_2 = linspace(0,2500e-6)

plot(d,N,'x')
plot(x_1,exp( m_1.nominal_value*x_1+b_1.nominal_value))
plot(x_2, exp(m_2.nominal_value*x_2+b_2.nominal_value))
yscale("log")
show()
