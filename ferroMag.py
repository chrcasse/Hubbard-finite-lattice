# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 08:50:04 2018

@author: Ch

Calculates the individual magnetization of each lattice point
by fix-point-iteration. 
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

ind = diag.indices(int(NNN**0.5)) #set of indices belonging to a certain sub-lattice


#Calculates the individual magnetization (mag), electron density (nn) and energy
def mag(m, U, n, nnn):
    ma = [0 for i in range(NNN)]
    nn = [0 for i in range(NNN)]
    en=0
    i=0
    x = linalg.eig(np.array(diag.init_hamilton(NNN,nnn,U,m)))
    ew1 = x[0].real
    ev1 = x[1]
    state1 = diag.good_order(ew1, ev1)
    y = linalg.eig(np.array(diag.init_hamilton(NNN,nnn,U,[-z for z in m])))
    ew2 = y[0].real
    ev2 = y[1]
    state2 = diag.good_order(ew2, ev2)
    for d in state1:
        nn = [nn[x]+np.abs(d[1][x])**2 for x in range(NNN)]
        if sum(nn) + 1e-3 > n*NN:
            break
        ma = [ma[x]+np.abs(d[1][x])**2 for x in range(NNN)]
        en += d[0]
        while i < len(state2) and d[0] > state2[i][0] and sum(nn) + 1e-3 < n*NN:
            ma = [ma[x]-np.abs(state2[i][1][x])**2 for x in range(NNN)]
            nn = [nn[x]+np.abs(state2[i][1][x])**2 for x in range(NNN)]
            en += state2[i][0]
            i+=1
    while i < len(state2) and sum(nn) + 1e-3 < n*NN:
        ma = [ma[x]-np.abs(state2[i][1][x])**2 for x in range(NNN)]
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
def ferro_mag(U, n,m, nn):
    mm = [1.2 for i in range(NNN)]
    en = 0
    en1 = 0
    z=0
    y=0
    while  or abs(en-en1) > 0.01:
        z += 1
        if z > 30:
#            print(abs(en-en1))
            if sum([(m[x]-mm[x])**2 for x in range(NNN)])/NN**2 < 0.0001 and abs(en-en1) > 0.01:
                m = mixing2(m, [ra.gauss(n_right(n), 0.1) for i in range(NNN)])
            y += 1
            z=0
            if y > 3:
                print(">>{}".format(U))
                break
        mm = m
        
        x = mag(mm, U, n, nn)
        m=x[0]
        nn=x[1]
        en = x[2]
        x = mag(m, U, n, nn)
        m=x[0]
        nn=x[1]
        en1 = x[2]
    
    return [m, nn, en1]

#functions for mixing of sets
def mixing2(a,b):
    return [0.7*a[i]+0.3*b[i] for i in range(len(a))]

def mixing1(a,b,c,d):
    return [0.5*a[i]+0.3*b[i]+0.1*c[i]+0.1*d[i] for i in range(len(a))]

def mixing(a,b,c,d):
    return [0.3*a[i]+0.4*b[i]+0.3*c[i]+0.0*d[i] for i in range(len(a))]


def n_right(n):
    if n > 1:
        return 2-n
    return n
#main function: return ind. magnetization, electron density and energy for 
#on-site-inter. 0<U<12 and given n
#plots heat-map for U in U_den
def m_U(n, U_den):
    m = [[],[]]
    nn = []
    en = []
    for q in range(3):
        m[0].append(12.1)
        x = ferro_mag(12.1, n,[n_right(n)*(1)**(i not in ind) for i in range(NNN)], [n for i in range(NNN)])
        m[1].append(x[0])
        nn.append(x[1])
        en.append(x[2])
    for U in [x/10000 for x in range(120000,0,-1000)]: 
        m[0].append(U)
        x = ferro_mag(U, n, mixing([n_right(n)*(1)**(i not in ind) for i in range(NNN)], m[1][-1], m[1][-2],m[1][-3]), mixing([n for i in range(NNN)], nn[-1], nn[-2],nn[-3]))
        x = ferro_mag(U, n, mixing1(x[0], m[1][-1], m[1][-2],m[1][-3]), mixing(x[1], nn[-1], nn[-2],nn[-3]))
#        x = ferro_mag(U, n, m[1][-1],nn[-1])
        m[1].append(x[0])
        nn.append(x[1])
        en.append(x[2])
        if U in U_den:
            magn = [[m[1][-1][x*int(NNN**0.5)+y] for x in range(int(NNN**0.5))] for y in range(int(NNN**0.5))]
            ax = sns.heatmap(magn, linewidth=0.0,cmap ="coolwarm")   
            plt.show()
    return [m, nn, en]

N = 1.0
print(N)
#
for N in [0.9,1.0,1.1]:#[x/20 for x in range(20,21)]:
    print(N)
    m = m_U(N,[9.0,2.0])
#    """m(U)&Energy"""
#    fobj = open("C:\\Users\\Ch\\Desktop\\uni\\6. Semester\\BA Störstelle\\Ferro144\\f_"+str(N)+".dat", "w")
#    fobj_en = open("C:\\Users\\Ch\\Desktop\\uni\\6. Semester\\BA Störstelle\\Ferro144\\energie_f_"+str(N)+".dat", "w")
#    for i in range(len(m[0][0])):
#        fobj.write(str(m[0][0][i]) + "\t" + str(sum(m[0][1][i])/NN) + "\n")
##        fobj.write(str(m[0][0][i]) + "\t" + str(sum([m[0][1][i][x] for x in ind])/NN1) + "\t" + str(sum([m[0][1][i][x] for x in list(set(range(NNN))-set(ind))])/NN2) + "\n")
#        fobj_en.write(str(m[0][0][i]) + "\t" + str(m[2][i]) + "\n")
##        print(sum(m[0][1][i])/NN)
##    print("\a")
#    fobj.close()
#    fobj_en.close()
#    """Standardabweichung von m(U)"""
#    fobj = open("C:\\Users\\Ch\\Desktop\\uni\\6. Semester\\BA Störstelle\\Ferro\\fs_stddev"+str(N)+".dat", "w")
#    for i in range(len(m[0][0])):
#        fobj.write(str(m[0][0][i]) + "\t" + str( np.std(m[0][1][i][0:55]+m[0][1][i][55:NNN])) + "\n")
#    fobj.close()
#    plt.plot(m[0][0] , [sum(m[0][1][i])/NN for i in range (len(m[0][0]))])
#    plt.show()
#    """Histogramm m(U)"""
#    fobj = open("hist"+str(N)+".dat", "w")
#    hist = np.histogram(m[0][1][50], bins=[y/100 for y in range(-100, 100)])[0]
#    j=0
#    for i in [y/100 for y in range(-100, 99)]:
#        fobj.write("{} {}\n".format(i, hist[j]))
#        j += 1
#    fobj.close()
    plt.plot(m[0][0] , m[2])
    plt.show()
print("\a")





x = ferro_mag(9,1.0,[1 for i in range(NNN)],[1 for i in range(NNN)])

    
    






