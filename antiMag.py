# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 08:50:04 2018

@author: Ch
"""
import diag  
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
import random as ra
import seaborn as sns

sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})

NNN=144 # Number of lattice points
NN=NNN-0 
NN1=NNN/2 # Number of sublattice points
NN2=NNN/2


ind = diag.indices(int(NNN**0.5))
dum_ind = ind

#Calculates the individual magnetization (mag), electron density (nn) and energy
def mag(m, U, n, nnn):
    #NN = 100
    ma = [0 for i in range(NNN)]
    nn = [0 for i in range(NNN)]
    en=0
    i=0
    x = linalg.eig(np.array(diag.init_hamilton_anti(NNN,nnn,U,m)))
    ew1 = x[0].real
    ev1 = x[1]
    state1 = diag.good_order(ew1, ev1)
    y = linalg.eig(np.array(diag.init_hamilton_anti(NNN,nnn,U,[-x for x in m])))
    ew1 = y[0].real
    ev1 = y[1]
    state2 = diag.good_order(ew1, ev1)
    #NN-=2
    for d in state1:
        nn = [nn[x]+np.abs(d[1][x])**2 for x in range(NNN)]
        if sum(nn)/NN  > n:
            break
        ma = [ma[x]+(-1)**(x not in ind)*np.abs(d[1][x])**2 for x in range(NNN)]
        en += d[0]
        while d[0] > state2[i][0] and sum(nn)/NN < n:
            ma = [ma[x]-(-1)**(x not in ind)*np.abs(state2[i][1][x])**2 for x in range(NNN)]
            nn = [nn[x]+np.abs(state2[i][1][x])**2 for x in range(NNN)]
            en += state2[i][0]
            i+=1
    while sum(nn)/NN  < n:
        ma = [ma[x]-(-1)**(x not in ind)*np.abs(state2[i][1][x])**2 for x in range(NNN)]
        nn = [nn[x]+np.abs(state2[i][1][x])**2 for x in range(NNN)]
        en += state2[i][0]
        i+=1
    for i in range(NNN):
        en -= U*(nn[i]-ma[i])*(nn[i]+ma[i])/4
    return [ma, nn, en/NN]

#Fix-point-iteration for on-site-interaction U,
#average electron density per site,
#individual magnetization m,
#individual electron density nn
#Returns individual magnetization, electron density and energy 
def anti_mag(U, n,m, nn,d):
    mm = [0.9 for i in range(NNN)]
    en1=0
    z=0
    y=1
    x = mag(m, U, n, nn)
    m=x[0]
    nn=x[1]
    en = x[2]
    while sum([(m[x]-mm[x])**2 for x in range(NNN)])/NN**2 > 0.0001:
        z+=1
        mm = m
        if z > 30:
            if y > 1:
                print(">>{}".format(U))
                break
#           
            y += 1
            z=0
            
        
        x = mag(m, U, n, nn)
        m=x[0]
        nn=x[1]
        en = x[2]

    return [m, nn, en]


def n_right(n):
    if n > 1:
        return 2-n
    return n

#functions for mixing of sets
def mixing3(a,b):
    return [0.7*a[i]+0.3*b[i] for i in range(len(a))]
    
def mixing(a,b,c,d):
    return [0.6*a[i]+0.3*b[i]+0.1*c[i]+0.0*d[i]+0.0*(1)**(i not in ind) for i in range(len(a))] # n > 1

def mixing2(a,b,c,d,e):
    return [0.4*a[i]+0.3*b[i]+0.3*c[i]+0.0*d[i]+0.0*e[i]+0.0*(1)**(i not in ind) for i in range(len(a))]

#main function: return ind. magnetization, electron density and energy for 
#on-site-inter. 0<U<12 and given n
#plots heat-map for U in U_den
def m_U(n, U_den):
    m = [[],[]]
    nn = []
    en = []
    ii = 0
    x = anti_mag(12.2, n,[1 for i in range(NNN)], [n for i in range(NNN)],0.0001)
    for q in range(4):
        #m[0].append(15.4-q/10)
        m[0].append(12.2)
        
        m[1].append(x[0])
        nn.append(x[1])
        en.append(x[2])
#        ii+=1
    for U in [x/10000 for x in range(120000,0,-1000)]:
#        print(U)
        m[0].append(U)
        x = []
        d=0.0002
        while x == []:
             x = anti_mag(12.2, n,[1 for i in range(NNN)], [n for i in range(NNN)],0.0001)
            #x  = anti_mag(U, n, mixing([n_right(n) for i in range(NNN)], m[1][-1], m[1][-2], m[1][-3]), mixing([n for i in range(NNN)], nn[-1],nn[-2],nn[-3]),d)
        
  
#            x  = anti_mag(U, n, mixing2(x[0], m[1][-1], m[1][-2], m[1][-3], m[1][-4]), mixing2(x[1], nn[-1],nn[-2],nn[-3],nn[-4]),d)
        
        m[1].append(x[0])
        nn.append(x[1])
        en.append(x[2])
        if U in U_den:
            magn = [[m[1][-1][x*int(NNN**0.5)+y]*(-1)**(x*int(NNN**0.5)+y not in ind) for x in range(int(NNN**0.5))] for y in range(int(NNN**0.5))]
            ax = sns.heatmap(magn, linewidth=0.0,cmap = "coolwarm")            
            plt.show()
#        print("{}:{}".format(U,sum(x[0])/NN))
            
        ii+=1
    return [m, nn, en]



for N in [0.9,1.0,1.1]:#[x/20 for x in range(21,21)]:#31:50
    print(N)
    m = m_U(N,[7,2])
    fobj = open("C:\\Users\\Ch\\Desktop\\uni\\6. Semester\\BA Störstelle\\Anti\\a144_"+str(N)+".dat", "w") # u wie gestoert
    fobj_en = open("C:\\Users\\Ch\\Desktop\\uni\\6. Semester\\BA Störstelle\\Anti\\energie_a144_"+str(N)+".dat", "w")
    for i in range(0, len(m[0][0])):
#        fobj.write(str(m[0][0][i]) + "\t" + str(sum(m[0][1][i])/NN) + "\n")
        fobj.write(str(m[0][0][i]) + "\t" + str(sum([m[0][1][i][x] for x in ind])/NN1) + "\t" + str(sum([m[0][1][i][x] for x in list(set(range(NNN))-set(ind))])/NN2) + "\n")
        fobj_en.write(str(m[0][0][i]) + "\t" + str(m[2][i]) + "\n")
#        print(sum(m[0][1][i])/NN)
#    print("\a")
    fobj.close()
    fobj_en.close()
    """Standardabweichung von m(U)"""
#    fobj = open("C:\\Users\\Ch\\Desktop\\uni\\6. Semester\\BA Störstelle\\Anti2\\fs_stddev"+str(N)+".dat", "w")
#    for i in range(len(m[0][0])):
#        fobj.write(str(m[0][0][i]) + "\t" + str( np.std([m[0][1][i][x] for x in dum_ind])) + "\n")
#    fobj.close()
    plt.plot(m[0][0] , [sum([m[0][1][i][x] for x in ind])/NN1 for i in range (len(m[0][0]))])
    plt.show()
print("\a")
    








