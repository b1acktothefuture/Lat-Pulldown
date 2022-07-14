import math
from fpylll import *
import numpy as np

def sign(x):
    return int(x>=0)

def norm(b):
    pass

def det(B):
    pass

def VolSphere(t):
    pass

def log2(x):
    return math.log(x)/math.log(2)

def ENUM(B, R, n):
    [c,u,y,v,sq] = [[0]*n for _ in range(5)]
    t = 0
    t_max = 0

    c1_bar = norm(B[0])  
    c[0] = c1_bar

    while(t<n):
        c[t] = c[t+1] + ((u[t] + y[t])**2)*(R[t][t]**2)
        if(c[t] < c1_bar and t>1):
            t = t - 1
            v[t] = 1
            y[t] = 0
            for i in range(t+1,t_max+1):
                y[t] += u[i]*R[t][i]/R[t][t]
            u[t] = -1*round(y[t])
            sq[t] = sign(u[t] - y[t])
        else:
            if(c[t] < c1_bar and t == 1):
                c1_bar = c1_bar
                b = 0
                pass
            else:
                pass
    return b


def NewENUM(B, R, s):
    sol =[]

    [n,m] = [len(B[0]),len(B)]
    L = []
    [t, t_max] = 1
    [c,u,y,sig, v] = [[0]*(n+1),[0]*(n+1),[0]*(n+1),[0]*(n+1),[0]*(n+1)]
    u[1] = 1
    c[1] = R[1][1]**2
    A = (n/4)*(det(B))**(1/n)
    while(t<= n):
        c[t] = c[t+1] + ((u[t] - y[t])*R[t][t])**2
        if(c[t]>=A):
            break
        # compute expected number of lattice points in Bt-1 using gaussian heuristic
        q_t = (A - c[t])**0.5
        volL = 1
        for i in range(1,t):
            volL *= R[i][i]
        beta_t = VolSphere(t-1)*(q_t**(t-1))/volL

        if( t == 1):
            b = [] #Linear combination of basis weighted by u
            norm_b = norm(b)**2
            if(norm_b**2 < A):
                sol.append(b)
                A = norm_b-1
                break
        if (beta_t >= 2**(-s + math.ceil(log2(log2(t))))):
            t -= 1
            y[t] = 0
            for i in range(t+1, t_max+1):
                y[t] += u[i]*R[t][i]/R[t][t]
            u[t] = -1*round(y[t])
            sig[t] = sign(u[t] - y[t])
            v[t] = 1
            continue
        L.append([])



    pass