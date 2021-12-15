import numpy as np
import math 
import random
import matplotlib.pyplot as plt
from scipy.stats import poisson
import pandas as pandas
from scipy.stats import uniform
import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad


### Radius in Micro Meter ###
rA549=4.8
rAG1522=6.7
rHela=8.35
rNB1RGB=7.4
rA1722=8.15
rV79=5.3
rCHO=6.36

### Ng is taken constant as Nbp is out of scope ###
pi=3.14
Ng=2500
### Different Values of Ns for in Micrometer ###
nsA549=(1.2)*(10**(-3))
nsAG1522=(4.2)*(10**(-4))
nsHela=(2.2)*(10**(-4))
nsNB1RGB=(3.2)*(10**(-4))
nsA1722=(2.4)*(10**(-4))
nsV79=(7.2)*(10**(-4))
nsCHO=(4.26)*(10**(-4))

### Calculate Ng for all of them ###

NgA549=((nsA549*16*100*7.5)/math.pi)*1000
NgAG1522=((nsAG1522*16*144*10.6)/math.pi)*1000
NgHela=((nsHela*16*219*13)/math.pi)*1000
NgNB1RGB=((nsNB1RGB*16*172*11.6)/math.pi)*1000
NgA1722=(nsA1722*16*209*12.7)/math.pi
NgV79=((nsV79*16*88*8.2)/math.pi)*1000
NgCHO=((nsCHO*16*127*9.9)/math.pi)*1000


### Nr and Lethality Criterion ###
N_r=np.random.uniform(0.04,0.08)
print(N_r)
Lambda=0.25 
### Probability of SSB converting to a Double Strand Break after ###


Pl_r =Lambda*(1-poisson.pmf(1,N_r)-poisson.pmf(0,N_r)-poisson.pmf(2,N_r))
sig=0.371

doseCHO = np.arange(0,6,0.0001)
doseCommon = np.arange(0,4,0.0001)

Se1= 100
Se1H=25
Se2=122
Se2H=17
Se3=54
Se3H=13
Se4=105
Se4H=78
Se5=120
Se5H=100
Se6=150
Se6H=72


X0=0.35
X1=0.04

Yl=((math.pi)*sig*Ng*doseCHO)/(16*Se6)
print(doseCHO)
print(Yl)

# plt.plot(Yl,dose)
# plt.show()
Pie=[]
Pie=[]
for i in Yl:
    if(i<(X0/X1)):
        Pie_surv=math.exp(-(1-X0)*i - X1*i*i)
        Pie.append(Pie_surv)
    else:
        Pie_surv=math.exp(-i)
        Pie.append(Pie_surv)

# print(Pie)
plt.plot(doseCHO,Pie)
plt.legend()
plt.show()



## Probability of Damage ###
damage=[]
for i in Pie:
    val=1-i
    damage.append(val)

plt.plot(doseCHO,damage)
plt.xlabel('dose')
plt.ylabel('Probability of damage')
plt.show()

### Listing The values of LET for Alln Cell Lines and also Hypoxic Conditions


def Yield  (Ng,Se):
    return (math.pi*sig*Ng*doseCHO)/(16*Se)

YlCHO=Yield(NgCHO,Se6)
print(max(YlCHO))

# YlCHO=(math.pi*sig*NgCHO)/(16*Se6)
# Ylhypoxic=(math.pi*sig*NgCHO)/(16*Se6H)

a=((1-X0)*math.pi*sig*NgCHO/(16*Se6))
b=X1*pow((math.pi*sig*NgCHO)/(16*Se6),2)

# ah=(1-X0)*Ylhypoxic
# bh=X1*pow(Ylhypoxic,2)

def func(n):
    return (math.exp(-a*n  -b*pow(n,2)))

y=list(map(func,doseCHO))

plt.yscale('log')
plt.plot(doseCHO,y,label='Cell survival deviation without repair')
plt.legend()
plt.show()

print(NgCHO)

