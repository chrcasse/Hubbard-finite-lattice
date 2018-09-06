# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 23:00:00 2018

@author: Ch
"""
import diag  
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
import random as ra
import seaborn as sns
#implementation of the fermidistribution
def fermi(T,mu,E):
    return (np.exp((E-mu)/T)+1)**-1

NNN=64
NN=NNN-0
NN1=NNN/2-0
NN2=NNN/2-0 
ind = diag.indices(int(NNN**0.5))
#entropy for a system of non-interacting fermions
def entropie(U,T,mu, energies):
    return -sum([fermi(T,mu,en)*np.log(fermi(T,mu,en))+(1-fermi(T,mu,en))*np.log(1-fermi(T,mu,en)) for en in energies])


def mag(m, U, n, nnn, T):
    ma = [0 for i in range(NNN)]
    nn = [0 for i in range(NNN)]
    en=0
    i=0
    x = linalg.eig(np.array(diag.init_hamilton_anti(NNN,nnn,U,m)))
    ew1 = x[0].real
    ev1 = x[1]
    state1 = diag.good_order(ew1, ev1)
    y = linalg.eig(np.array(diag.init_hamilton_anti(NNN,nnn,U,[-z for z in m])))
    ew2 = y[0].real
    ev2 = y[1]
    state2 = diag.good_order(ew2, ev2)
    mu = diag.good_order2(list(ew1)+list(ew2))[int(NN*n)]
    for d in state1:
        nn = [nn[x]+np.abs(d[1][x])**2*fermi(T, mu, d[0]) for x in range(NNN)]
        ma = [ma[x]+(-1)**(x not in ind)*np.abs(d[1][x])**2*fermi(T, mu, d[0]) for x in range(NNN)]
        en += d[0]*fermi(T, mu, d[0])
        
        
    for d in state2:
        nn = [nn[x]+np.abs(d[1][x])**2*fermi(T, mu, d[0]) for x in range(NNN)]
        ma = [ma[x]-(-1)**(x not in ind)*np.abs(d[1][x])**2*fermi(T, mu, d[0]) for x in range(NNN)]
        en += d[0]*fermi(T, mu, d[0])
    for i in range(NNN):
        en -= U*(nn[i]-ma[i])*(nn[i]+ma[i])/4
    
    return [ma, nn, (en-T*entropie(U,T,mu,list(ew1)+list(ew2)))/NN]

#Fix-point-iteration for on-site-interaction U,
#average electron density per site,
#individual magnetization m,
#individual electron density nn
#Returns individual magnetization, electron density and energy 
def anti_mag(U, n,m, nn,T):
    mm = [1.2 for i in range(NNN)]
    en = 0
    z=0
    while sum([(m[x]-mm[x])**2 for x in range(NNN)])/NN**2 > 0.00005:
        z += 1
        if z > 50:
            print(">>{}".format(sum(m)/NN))
#            m = [ra.gauss(n*(-1)**(i not in ind), 0.1) for i in range(NNN)]
            return [m, nn, 0]
            z=0
        mm = m
        
        x = mag(mm, U, n, nn,T)
        m=x[0]
        nn=x[1]
        en = x[2]
#    print(">>>{}".format(sum([(m[x]-mm[x])**2 for x in range(NNN)])/NN))#sum([(m[x]-mm[x])**2 for x in range(NNN)]))
        
        if sum(m) < 0:
            m = [ra.gauss(n, 0.1) for i in range(NNN)]
    return [m, nn, en]

def mixing1(a,b,c):
    return [abs(0.2*a[i]+0.5*b[i]*0.3*c[i])*(-1)**(i not in ind) for i in range(len(a))]

def mixing(a,b,c,d):
#    return [0.0*a[i]+0.5*b[i]+0.5*c[i]+0.0*d[i] for i in range(len(a))]
    return [0.0*a[i]+0.5*b[i]+0.5*c[i]+0.0*d[i] for i in range(len(a))] # anti bis 15
#ferro: 50/50

def n_right(n):
    if n > 1:
        return 2-n
    return n
#main function: return ind. magnetization, electron density and energy for 
#temerature 0<T<3, given U and given n
def m_T(n, U):
    m = [[],[]]
    nn = []
    en = []
    x = anti_mag(U, n,[n_right(n)*(1)**(i not in ind) for i in range(NNN)], [n for i in range(NNN)],1.1)
    for q in range(3):
        m[0].append(0.01)
        
        m[1].append(x[0])#[n_right(n)*(1)**(i not in ind) for i in range(NNN)])
        nn.append(x[1])
        en.append(x[2])
    for T in [x/100000 for x in range(1000, 300000, 100)]:
        print(T)
        m[0].append(T)
        x = anti_mag(U, n, mixing([n_right(n)*(1)**(i not in ind) for i in range(NNN)], m[1][-1], m[1][-2],m[1][-3]), mixing([n for i in range(NNN)], nn[-1], nn[-2],nn[-3]),T)
#        x = ferro_mag(U, n, m[1][-1],nn[-1])
        m[1].append(x[0])
        nn.append(x[1])
        en.append(x[2])
        
    return [m, nn, en]

N = 1.1
print(N)
#
for U in [10]:
    print(U)
    m = m_T(N,U)
#    fobj = open("C:\\Users\\Ch\\Desktop\\uni\\6. Semester\\BA Störstelle\\AntiT\\a8_"+str(N)+".dat", "w")
#    fobj_en = open("C:\\Users\\Ch\\Desktop\\uni\\6. Semester\\BA Störstelle\\AntiT\\energie_a8_"+str(N)+".dat", "w")
#    for i in range(len(m[0][0])):
#        fobj.write(str(m[0][0][i]) + "\t" + str(sum([m[0][1][i][x] for x in ind])/NN1) + "\n")
##        fobj.write(str(m[0][0][i]) + "\t" + str(sum([m[0][1][i][x] for x in ind])/NN1) + "\t" + str(sum([m[0][1][i][x] for x in list(set(range(NNN))-set(ind))])/NN2) + "\n")
#        fobj_en.write(str(m[0][0][i]) + "\t" + str(m[2][i]) + "\n")
##        print(sum(m[0][1][i])/NN)
##    print("\a")
#    fobj.close()
#    fobj_en.close()
    plt.plot(m[0][0] , [sum([m[0][1][i][x] for x in ind])/NN1 for i in range (len(m[0][0]))])
    plt.show()
print("\a")

    
    






