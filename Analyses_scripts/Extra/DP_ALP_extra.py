# ---------- Libaries ----------
#from matplotlib.lines import _LineStyle
import numpy as np
import pandas as pd
from numpy.ma.core import sqrt
import scipy.stats as stats
import math
import sympy
from scipy.special import binom

from sympy.parsing.mathematica import parse_mathematica
from sympy.parsing.mathematica import mathematica
import scipy.integrate as integrate
from sympy import var
from sympy import N

def f1(x,a_n,b_n,c_n,d_n,a_d,b_d,c_d):
    return ( a_n*x**3 + b_n*x**2 + c_n*x +d_n )/(a_d*x*(b_d*x+c_d)**2)

def pol_cos(x,a,b,c,d):
    return a*np.cos(b*x+d)+c

def pol_const(x,c):
    return c

def percentage_error(data,total):
    if total != 0 and data != 0:
        return (data/total)*sqrt(1/data + 1/total)
    else:
        return 0

def append_func_zero(listTOT,list1,list2):
    if float(list2.sum()) == 0:
        listTOT = np.append(listTOT,0.0)
    else:
        #print(np.divide(list1,list2).sum()/len(list1))
        #print(list1.sum()/list2.sum())
        listTOT = np.append(listTOT,list1.sum()/list2.sum())
    return listTOT

def for_back_calcs(res_for,res_back,errs_for,errs_back,indi,tots,cos_values,cos_min,cos_max):
    # --- Forwards ---
    temp1 = np.array([ indi[c_val] for c_val in range(0,len(cos_values)) if cos_values[c_val] >= 0.0 and cos_values[c_val] <= cos_max ])
    temp2 = np.array([ tots[c_val] for c_val in range(0,len(cos_values)) if cos_values[c_val] >= 0.0 and cos_values[c_val] <= cos_max ])
    res_for = append_func_zero(res_for,temp1,temp2)
    errs_for = np.append(errs_for,percentage_error(temp1.sum(),temp2.sum()))

    # --- Backwards ---
    temp1 = np.array([ indi[c_val] for c_val in range(0,len(cos_values)) if cos_values[c_val] >= cos_min and cos_values[c_val] <= 0.0 ])
    temp2 = np.array([ tots[c_val] for c_val in range(0,len(cos_values)) if cos_values[c_val] >= cos_min and cos_values[c_val] <= 0.0 ])
    res_back = append_func_zero(res_back,temp1,temp2)
    errs_back = np.append(errs_back,percentage_error(temp1.sum(),temp2.sum()))

    return res_for, res_back, errs_for, errs_back

def sigmod_func(x,k):
    if (1+abs(x)**k)**(float(1)/float(k)) != 0:
        return x/((1+abs(x)**k)**(float(1)/float(k)))
    else:
        return 0

def R_func(x_vals,mA,Eng1,Eng2,me,line_str,s_denom,s_num,s_other):
    out = np.array([])
    gE,gA,gX,Eng1,Eng2,me= var('gE gA gX Eng1 Eng2 me')

    Eng1=4
    Eng2=7
    gA=1
    gE=10**4
    gX=1
    me=0.511*10**(-3)

    for x in x_vals:       
        if len(s_other) != 0 and pd.isna(eval(s_other)):
            out_i = 1.0
        else:
            out_denom = eval(s_denom) 
            out_num = eval(s_num)

            if math.isnan(out_denom):
                out_i = 0.0
            elif math.isnan(out_num):
                out_i = 1.0
            elif np.round(out_denom,5) == 0.0 or np.round(out_denom,5) == -0.0 or np.round(out_denom,5) == np.round(out_num,5):
                if x == -1.0:
                    x = -0.999
                elif x == 1.0:
                    x = 0.999
                out_i = eval(line_str)
            else:
                out_i = eval(line_str)
        
        #out_i = eval(line_str)
        out = np.append(out,out_i)

    return out

def split_fraction(line_ratio):
    line_str = str(line_ratio).replace("sign","np.sign").replace("Abs","np.abs").replace("sin","np.sin").replace("cos","np.cos").replace("anp.cos","np.arccos").replace("sqrt","np.sqrt")
    
    s_denom = line_str.split("/")[-1]  # denominator
    out_denom = 0
        
    s_num = "" # numerator
    out_num = 0

    s_other = "" # other fractions
    out_other = 0
    
    for i in range(0,len(str(line_ratio).split("/"))):
        if i == 0:
            s_num = line_str.split("/")[i]
        elif i != (len(str(line_ratio).split("/"))-1):
            s_num = s_num + "/" + line_str.split("/")[i]

        if i != 0 and i != (len(str(line_ratio).split("/"))-1):
            test1=line_str.split("/")[i]
            test2=test1.split(") ")
            test3=test2[0]+")"

            if str(test3).count("(") < str(test3).count(")"): # uneven parentises
                n_test = str(test3).count(")") - str(test3).count("(")
                test4=str(test3)[:-n_test]
                s_other = test4
            else:
                s_other = str(test3)
    return line_str, s_denom, s_num, s_other