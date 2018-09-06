# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 18:34:41 2018

@author: Ch
"""

import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt


def get_indexI(y, nn, i):
    if i % int(np.sqrt(nn))== 0:
        return i - int(np.sqrt(nn))
    
    return i

def get_indexII(y, nn, i):
    if i % int(np.sqrt(nn)) == int(np.sqrt(nn))-1:
        return i + int(np.sqrt(nn))
    return i

def get_indexIII(nn, i):
    if i >= nn:
        return i - nn 
    
    return i

def get_indexIV(nn, i):
     
    if i < 0:
        return i + nn
    return i
    

def init_hamilton_anti(nn, ne, U, mm):
    h = []
    ind = indices(int(nn**0.5))
    for i in range(nn):
        h.append([0 for j in range(nn)])
    for y in range(nn):
        h[y][y] = U*(ne[y]-(-1)**(y not in ind)*mm[y])/2
        h[y][get_indexI(y, nn, y+1)] -= 1
        h[y][get_indexII(y, nn, y-1)] -= 1
        
        h[y][get_indexIII(nn, y+int(nn**0.5))] -= 1
        h[y][get_indexIV(nn, y-int(nn**0.5))] -= 1
        
        
#    for s in [0]:
#        h[s] = [0 for j in range(nn)]
#        for y in range(nn):
#            h[y][s] = 0
#        h[s][s] = U*(ne[s]-mm[s])/2

    
    return h

def init_hamilton_anti_einfach(nn, ne, U, mm):
    h = []
    ind = indices(int(nn**0.5))
    for i in range(nn):
        h.append([0 for j in range(nn)])
    for y in range(nn):
        h[y][y] = U*(ne-(-1)**(y not in ind)*mm)/2
        h[y][get_indexI(y, nn, y+1)] -= 1
        h[y][get_indexII(y, nn, y-1)] -= 1
        
        h[y][get_indexIII(nn, y+int(nn**0.5))] -= 1
        h[y][get_indexIV(nn, y-int(nn**0.5))] -= 1
        
#        
    for s in [55]:
        h[s] = [0 for j in range(nn)]
        for y in range(nn):
            h[y][s] = 0
        h[s][s] = U*(ne-(-1)**(s not in ind)*mm)/2

    
    return h

def init_hamilton(nn, ne=1, U=0, mm=0):
    h = []
    #nn+=2
    for i in range(nn):
        h.append([0 for j in range(nn)])
    for y in range(nn):
        h[y][y] = U*(ne[y]-mm[y])/2
#        h[y][y] = U*(ne-mm)/2
        h[y][get_indexI(y, nn, y+1)] += -1
        h[y][get_indexII(y, nn, y-1)] += -1
        
        h[y][get_indexIII(nn, y+int(nn**0.5))] += -1
        h[y][get_indexIV(nn, y-int(nn**0.5))] += -1
#        
    for s in [55]:
        h[s] = [0 for j in range(nn)]
        for y in range(nn):
            h[y][s] = 0
        h[s][s] = U*(ne[s]-mm[s])/2

    return h

def init_hamilton_einfach(nn, ne=1, U=0, mm=0):
    h = []
    #nn+=2
    for i in range(nn):
        h.append([0 for j in range(nn)])
    for y in range(nn):
        h[y][y] = U*(ne-mm)/2
        h[y][get_indexI(y, nn, y+1)] += -1
        h[y][get_indexII(y, nn, y-1)] += -1
        
        h[y][get_indexIII(nn, y+int(nn**0.5))] += -1
        h[y][get_indexIV(nn, y-int(nn**0.5))] += -1
        
#    for s in [0,20,60,85]:
#        h[s] = [0 for j in range(nn)]
#        for y in range(nn):
#            h[y][s] = 0
#        h[s][s] = U*(ne-(1)**(s not in ind)*mm)/2
    return h

def init_hamilton_einfach_para(nn, ne=1, U=0, mm=0):
    h = []
    #nn+=2
    for i in range(nn):
        h.append([0 for j in range(nn)])
    for y in range(nn):
        h[y][y] = U*(ne)/2
        h[y][get_indexI(y, nn, y+1)] += -1
        h[y][get_indexII(y, nn, y-1)] += -1
        
        h[y][get_indexIII(nn, y+int(nn**0.5))] += -1
        h[y][get_indexIV(nn, y-int(nn**0.5))] += -1
        
    for s in [55]:
        h[s] = [0 for j in range(nn)]
        for y in range(nn):
            h[y][s] = 0
        h[s][s] = 1e10#U*(ne-(-1)**(s not in ind)*mm)/2

    return h



def indices(N):
    dum = []
    for x in range(0,N,2):
        for i in range(x*N,x*N+N,2):
            dum.append(i)
        for i in range(x*N+N+1,x*N+2*N,2):
            if i <= N**2:
                dum.append(i)
    return dum

#print(indices(10))

def norm(v):
    sq_sum = sum([np.abs(x)**2 for x in v])**0.5
    v = [x/sq_sum for x in v]
    return v

def get_key(item):
    return item[0]

def get_key2(item):
    return item

def good_order(ew, ev):
    dum = []
    for i in range(len(ew)):
        dum.append([ew[i], ev[:,i]])
    return sorted(dum, key=get_key)

def good_order2(ew):
    dum = []
    for i in range(len(ew)):
        dum.append(ew[i])
    return sorted(dum, key=get_key2)





#for N in [2500]:
#    ind = indices(N)
#    print(N)
#    x = linalg.eig(np.array(init_hamilton_einfach(N,1,12,1))) #anti:0.8-12
#    ew = x[0]
#    ev = x[1]
#    states = good_order(ew, ev)
#    i=0
##    for d in states:
##        print("{}>>>{}".format(i,d[0].real))
##        i+=1
#    dos = [[],[]]
#    i = 0
#    fobj = open("dos_f"+str(N)+".dat", "w")
#    hist = np.histogram(ew, bins=[y/6 for y in range(-40, 260)], density=True)[0]
#    plt.plot([y/6 for y in range(-40, 259)],hist)
#    plt.show()
#    j=0
#    z = [round(s,3).real for s in ew]
#    
#    ew = list(set(z))
#    for i in ew:
#        
#        fobj.write("{} {}\n".format(i.real, z.count(i.real)/N))
#        j += 1
#    fobj.close()




