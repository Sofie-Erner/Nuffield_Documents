import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors
plt.style.use('tableau-colorblind10')
import numpy as np
from pathlib import Path
import sys
import os
import math
import statistics

from sklearn.preprocessing import StandardScaler 
from sklearn.preprocessing import MinMaxScaler

from scipy.stats import poisson

def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)

    if b > 1:
        if a == "1":
            return r'$10^{{{}}}$'.format(b)
        else:
            return r'${} \times 10^{{{}}}$'.format(a,b)
    elif b == -1:
        return r'0'
    elif b < 0:
        return r'$10^{{{}}}$'.format(b)
    else:
        return r'{}'.format(round(x))

def plot_min_max(z_vals):
    maxA, maxB = '{:.0e}'.format(np.max(z_vals)).split('e')
    setMax = np.max(z_vals)

    if int(maxB) > 0:
        setMin = 0.1
    elif int(maxB) == -1:
        setMin = 0.01
        setMax = 1
    else:
        setMin = 10**(int(maxB)-2)

    return [setMin, setMax]

def R_func(x_vals,str_fit):
    out = np.array([])
    for x in x_vals:       
        out_i = eval(str_fit)

        out = np.append(out,out_i)
    return out

def corr_precision(val):
    expo = round(np.log10(val))

    if expo < -3:
        return 1*10**(-3)
    else:
        return 1*10**expo
    

def find_signal(N_SM,N_events,p_crit):
    i = 0
    ptemp = poisson.cdf(k=N_SM, mu=(N_SM + N_events))

    if round(ptemp,2) == p_crit:
        return i
    elif ptemp < p_crit:
        pVal = 1-ptemp
        cVal = 1-p_crit
    else:
        pVal = ptemp
        cVal = p_crit

        
    while round(pVal,2) > cVal:
        ptemp = poisson.cdf(k=N_SM, mu=(N_SM + (1+i)*N_events))
        #print(ptemp)
        if ptemp == 0:
            i-= 10**(-3)
            pVal = 1 - ptemp
            cVal = 1 - p_crit
        elif ptemp < 0.1:
            if round(np.log10(N_events)) > round(np.log10(N_SM)):
                i_temp = corr_precision(ptemp)*10**(-round(np.log10(N_SM)))

                while round(np.log10(i_temp)) < 3:
                    i_temp = i_temp*10
    
                i-= i_temp
            else:
                i-=corr_precision(ptemp)*10**(-round(np.log10(N_events)))
            pVal = 1 - ptemp
            cVal = 1 - p_crit
        else:
            i+=corr_precision(ptemp)*10**(-round(np.log10(N_events)))
            pVal = ptemp
            cVal = p_crit
    
    return i