from numpy import genfromtxt, array,linspace, exp, log,sqrt
from uncertainties import  ufloat
from matplotlib.pyplot import *

from Tools import lin_reg , make_LaTeX_table
from sympy.physics.units import length
data = genfromtxt("beta_data1.csv", delimiter=',', unpack="true")

N_3 = data[1]
d_3 = data[0]*1e-6*2700
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


print " m_1: %s b_1:%s" % (m_1,b_1)
print " m_2: %s b_2:%s" % (m_2,b_2)
x_1 = linspace(0,2.8)
x_2 = linspace(0,7)

R_max = (b_2-b_1)/(m_1-m_2)

print "R_Max: %s m " % (R_max)

R_max*= 0.1
E_Max = 1.92*(R_max**2+0.22*R_max)**0.5 
R_max*= 10
print "E_max : %s MeV" % E_Max
xlabel(r"Reichweite [$\frac{kg}{m^2}$]")
ylabel("$A$ in [$s^{-1}$]")
ylim(1,1000)

errorbar(d_3, A_3, A_3_err,0,'.', label= "Messwerte")
plot(x_1,exp( m_1.n*x_1+b_1.n), label= "Lineare Regession Teil A")
plot(x_2, exp(m_2.n*x_2+b_2.n), label= "Lineare Regession Teil B")

plot([R_max.n], [exp(m_1.n*R_max.n +b_1.n)],'ro', )
yscale("log")

legend()

savefig("beta_log.png")


# Tabelle erstellen

data = [[int(d_3[i]*1e6),ufloat(N_3[i],sqrt(N_3[i])), ufloat(A_3[i], A_3_err[i])] for i in range(len(d_3))]

header = [ r"$\frac{d}{\si{\micro\meter}}$",r"$N$" , r"$\frac{A}{\si{\second^{-1}}}$"]

#print make_LaTeX_table(array(data), header)

print (E_Max-0.314)/0.314