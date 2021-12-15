import numpy as np
import matplotlib.pyplot as plt
from random import randint
import math
import scipy. stats as ss
from math import exp


dose = np.arange(0,6,0.001)
x0 = 0.35
x1 = 0.04
sig = 0.371
#Ng = 2500
Se = 150
lamda=0.15
#Chirayu



Dn_A549= 9.6
Dn_AG1522= 13.4
Dn_Hela=16.7
Dn_NB1RGB=14.8
Dn_A1722=16.3
Dn_V79=10.6
Dn_CHO = 12.7

ns_A549=(1.2)*(10**(-3))
ns_AG1522=(4.2)*(10**(-4))
ns_Hela=(2.2)*(10**(-4))
ns_NB1RGB=(3.2)*(10**(-4))
ns_A1722=(2.4)*(10**(-4))
ns_V79=(7.2)*(10**(-4))
ns_CHO=(4.26)*(10**(-4))

Se1= 100
Se1H=25
Se2=122
Se2H=17
Se3=70
Se3H=15
Se4=54
Se4H=13
Se5=105
Se5H=78
Se6=120
Se6H=100
Se7=150
Se7H=72

#gradNumOfLesion = Ns*sig

def getZ(Dn):
    return (math.pi*Dn)/4
    
z_A549= getZ(Dn_A549)
z_AG1522= getZ(Dn_AG1522)
z_Hela=getZ(Dn_Hela)
z_NB1RGB=getZ(Dn_NB1RGB)
z_A1722=getZ(Dn_A1722)
z_V79=getZ(Dn_V79)
z_CHO = getZ(Dn_CHO)

def getNbp(Ns,An,z):
    return (160*An*z*Ns)/(3*math.pi)

Nbp_A549= getNbp(ns_A549,144,z_A549)
Nbp_AG1522= getNbp(ns_AG1522,100,z_AG1522)
Nbp_Hela=getNbp(ns_Hela,219,z_Hela)
Nbp_NB1RGB=getNbp(ns_NB1RGB,172,z_NB1RGB)
Nbp_A1722=getNbp(ns_A1722,209,z_A1722)
Nbp_V79=getNbp(ns_V79,88,z_V79)
Nbp_CHO = getNbp(ns_CHO,127,z_CHO)


def getNg(Nbp):
    return Nbp/3.33

Ng_A549= getNg(Nbp_A549)
Ng_AG1522= getNg(Nbp_AG1522)
Ng_Hela=getNg(Nbp_Hela)
Ng_NB1RGB=getNg(Nbp_NB1RGB)
Ng_A1722=getNg(Nbp_A1722)
Ng_V79=getNg(Nbp_V79)
Ng_CHO = getNg(Nbp_CHO)

def Yield  (Ng,Se):
    return (math.pi*sig*Ng*dose)/(16*Se)

YlA549 = Yield(Ng_A549,Se1)
YlAG1522 = Yield(Ng_AG1522,Se2)
YlHela = Yield(Ng_Hela,Se3)
YlNB1RGB = Yield(Ng_NB1RGB,Se4)
YlA1722 = Yield(Ng_A1722,Se5)
YlV79 = Yield(Ng_V79,Se6)
YlCHO = Yield(Ng_CHO,Se7)
 

# N_r=randint(1,10)
N_r=np.random.uniform(0.04,0.08)
print(N_r)


Pl_r= (1- (ss.poisson.pmf(0, N_r) + ss.poisson.pmf(1, N_r) + ss.poisson.pmf(2, N_r)))
Pl_r=Pl_r*lamda
print(Pl_r)

#actual cell survival curve
cell_survival=[]

for i in YlA549:
    
    if(i<(x0/x1)):
        cell_surv=exp(-(1-x0)*i - x1*i*i)
        cell_survival.append(cell_surv)
    else:
        cell_surv=exp(-i)
        cell_survival.append(cell_surv)

damage=[]
Pis=[]
for i in cell_survival:
    val=1-i
    damage.append(val)



plt.xlabel('dose')
plt.ylabel('Probability of damage')

#expected cell survival curve if cell does not have repair factor for small amount of dose
cell_surv1=[]
for i in YlA549:
    cell_surv=exp(-i)
    cell_surv1.append(cell_surv)
#expected cell survival curve if cell does  have repair factor for small amount of dose    

cell_surv2=[]
for i in YlA549:
    cell_surv=exp(-(1-x0)*i - x1*i*i)
    cell_surv2.append(cell_surv)







plt.plot(dose,cell_survival,color ='r',label ='actual survival curve')
plt.plot(dose,cell_surv1,color ='g',label ='expected survival curve in absence of repair factor of')
plt.xlabel('dose')
plt.ylabel('Probability of cell survival ')
plt.legend()
plt.show()

plt.plot(dose,cell_survival,color ='r',label ='actual survival curve')
plt.plot(dose,cell_surv2,color ='b',label ='expected survival curve due to repair factor of cell')
plt.xlabel('dose')
plt.ylabel('Probability of cell survival ')
plt.legend()
plt.show()


plt.xlabel('dose')
plt.ylabel('Probability of damage')
plt.plot(dose,damage)
plt.show()