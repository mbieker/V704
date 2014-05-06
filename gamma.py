from numpy import genfromtxt, sqrt, array, log, linspace, pi, exp
from uncertainties import ufloat
import uncertainties.umath as umath
from scipy.constants import *
from matplotlib.pyplot import *
data = genfromtxt("gamma_data1.csv", delimiter=',', unpack="true")

N_0 = 165.0
d_1 = data[0]*1e-3
N_1 = data[1]
d_2 = data[2]*1e-3
N_2 = data[3]

def fehler(N,N_0,t,t_0):
    return sqrt(N/t**2 + N_0/t_0**2)



def mu_theo(Z, M, rho):
    r_e = 2.82e-15
    E_gamma = 0.662e6 * electron_volt
    eps = E_gamma/(c**2*m_e)
    bracket1 = 2*(1.0+eps)/(1.0+2*eps)- log(1+2*eps)/eps
    bracket2 = ((1+eps)/eps**2)*bracket1 + 0.5*log(1+2*eps)/eps - (1+3*eps)/(1+2*eps)**2
    sigma = 2*pi*r_e**2*bracket2
    return (Z*N_A*rho*sigma)/M
t_0 = 900.0 
t_1 = 300.0
t_2 = 100.0

A_0 = N_0/t_0
# Messung fuer Kupfer

A_1 = N_1/t_1
A_1_err = array([fehler(i,N_0,t_1,t_0) for i in N_1])


#Daten Plotten

xlabel("Schichtdicke [$m$]")
ylabel("$A-A_0$ in [$s^{-1}$]")
xlim(0,0.02)
ylim(10,32)
errorbar(d_1,A_1-A_0,A_1_err,0, "x", label="Messwerte")
legend()
savefig("gamma1_lin.png")


#Lineare  Regression

from Tools import lin_reg, make_LaTeX_table

m_1, b_1 = lin_reg(d_1,log(A_1-A_0))


print "mu %s und A(0) = %s" % (m_1,b_1)
print "A_0:%s " % umath.exp(b_1)

#Plot LinFit

yscale("log")
ylim(10,100)
d = linspace(0,0.02)
plot(d,exp(m_1.n * d+b_1.n), label = "Lineare Regression")
legend()
savefig("gamma1_log.png")

#Latex Tabelle erzeugen

data = array([[int(d_1[i]*1e3),int(N_1[i]), ufloat(A_1[i]-A_0, A_1_err[i]) ]for i in range(7)])
#print make_LaTeX_table(data,  [r'{$\frac{d}{\si{\milli\meter}}$} ',r'$N$' ,r'{ $\frac{A-A_0}{\si{\second^{-1}}}$ }'])



#Theorie vgl.

mu_th  = mu_theo(29, 63.546, 8.92e3)
print "mu_theo_cu: %s " % mu_th







# Das Gleiche nochmal fuer Pb

A_2 = N_2/t_1
A_2_err = array([fehler(i,N_0,t_1,t_0) for i in N_1])


#Daten Plotten
close()
yscale("linear")
xlabel("Schichtdicke [$m$]")
ylabel("$A-A_0$ in [$s^{-1}$]")

errorbar(d_2,A_2-A_0,A_2_err,0, "x", label ="Messwerte")
legend()
savefig("gamma2_lin.png")


#Lineare  Regression


m_1, b_1 = lin_reg(d_2,log(A_2-A_0))


print "mu_2 %s und lnA_2(0) = %s" % (m_1,b_1)
print "A_0:%s " % umath.exp(b_1)

#Plot LinFit

yscale("log" )

d = linspace(0, 0.035)
plot(d,exp(m_1.n * d+b_1.n),label = "Lineare Regression")
xlim(0,0.035)
legend()
show()

#Latex Tabelle erzeugen

data = array([[int(d_2[i]*1e3),int(N_2[i]), ufloat(A_2[i]-A_0, A_2_err[i]) ]for i in range(7)])
#print make_LaTeX_table(data,  [r'{$\frac{d}{\si{\milli\meter}}$} ',r'$N$' ,r'{ $\frac{A-A_0}{\si{\second^{-1}}}$ }'])



#Theorie vgl.

mu_th  = mu_theo(82, 207.2, 11.342e3)
print "mu_theo_pb: %s " % mu_th

