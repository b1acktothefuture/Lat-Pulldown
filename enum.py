from fpylll import *
import numpy as np

def sign(x):
    return int(x>=0)

def norm(b):
    pass

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


    pass