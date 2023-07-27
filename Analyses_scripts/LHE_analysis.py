# ---------- Libaries ----------
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('tableau-colorblind10')

import sys
import os
import math
import numpy as np
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.special import expit
import re

from sympy import N

from Extra.stat_extra import *
from Extra.fit_funcs import *
from Extra.DP_ALP_extra import *

# ---------- Command Line Arguments ----------
if ( len(sys.argv) < 3 ): # less than three command line arguments
    print("Not correct number of command line arguments")
    quit()

# ---------- Variables ----------
direc = sys.argv[1] # directory path
type = sys.argv[2] # type of analysis
prob_text = "Percentage (%)" #"Probability"

if type == "contour":
    # ---------- SM values
    SM_tot_prob = 0.521
    SM_tot_err = 0.013018978653044034
    SM_for = 0.511
    SM_for_err = 0.019
    SM_back = 0.528
    SM_back_err = 0.018

    # ---------- Mass and coupling variables
    m_values = np.array([])
    g_values = np.array([])

    g_text = "$g$"
    m_text = "$m$"

    if ( len(sys.argv) < 4 ): # less than three command line arguments
        print("Not correct number of command line arguments")
        quit()
    else:
        proc = sys.argv[3] # type of process

        if proc == "DP":
            g_text = "$g_{A'}$"
            g_const = "0.001"
            m_text = "$M_{A'}$"
        elif proc == "ALP":
            g_text = "$f_{a}$"
            g_const = "40000.0"
            m_text = "$M_{ALP}$"
        
        m_const = "1.0"
elif type == "SM":
    SM_more = 0 # more than one SM type
    SM_proc = ""



data_no_cuts = {} # directory for data_no_cuts from histograms
data_w_cuts = {} # directory for data_no_cuts from histograms
ending='.dat'
titles = np.array([])  # empty array for titles of histograms
hist_names = np.array(["hsA","hsE","hsE1","hsEm1"])
n_h = len(hist_names)

cos_tots_w_cuts = {}
cos_Se_Sa_w_cuts = {}
cos_probs_w_cuts = {}

cos_tots_no_cuts = {}
cos_Se_Sa_no_cuts = {}
cos_probs_no_cuts = {}

cos_all_s = {}
cos_h_minus = {}

cos_vals = np.array([])
cos_min = np.cos(163*(np.pi/180))
cos_max = np.cos(30*(np.pi/180))

pattern = r'[0-9]'

# ---------- Read through files ----------
# return all files as a list
for file in os.listdir(direc):
    if file.endswith(ending) and "out" not in file:
        #print(file)
        in_file = open(os.path.join(direc, file),'r')
        lines = in_file.readlines()

        title = file.replace(ending,'')
        if type == "E":
            title = int(title.replace('_7', ''))
            titles = np.append(titles,title)

        elif type == "SM":
            title = title.replace('SM_','').replace('run_','')
            title = re.sub(pattern, '', title) # remove numbers
            a_title = title.split('_')
            if a_title[0] == "all":
                title = title
            elif a_title[0] != "ee":
                if "aa_aaa" in title:
                    title = title.replace("aa_aaa","$ \gamma \gamma (\gamma)$")
            else:
                title_temp = ""
                for i in range(1,len(a_title)):
                    title_temp = title_temp + a_title[i]
                title_temp = title_temp.replace('a','\gamma ').replace('nn','\\nu \\bar{\\nu} ').replace('ee','e^+ e^- ').replace('z','Z \\rightarrow')
                title = "$ "+ title_temp + " $"
            
            if SM_proc == "" and a_title[0] != "all":
                if "ee" in a_title[1]:
                    SM_proc = "e"
                elif "n" in a_title[1]:
                    SM_proc = "n"
                else:
                    SM_proc = "a"
            elif a_title[0] != "all":
                if "ee" in a_title[1]:
                    if SM_proc != "e":
                        SM_more = 1
                elif "n" in a_title[1]:
                    if SM_proc != "n":
                        SM_more = 1
                else:
                    if SM_proc != "a":
                        SM_more = 1

            title = title.replace('_','')
            titles = np.append(titles,title)

        elif type == "contour":
            a1_title = title.split("_") # split by _ to separate mA' and g1'
            m = float(a1_title[0].replace("m",""))
            g = float(a1_title[1].replace("g",""))

            if m_values.size == 0:
                m_values = np.array([m])
                g_values = np.array([g])
            else:
                if m not in m_values: # mA'
                    m_values = np.append(m_values,m)
                if g not in g_values: # g1'
                    g_values = np.append(g_values,g)
            title = str(m)+"_"+str(g)
            titles = np.append(titles,title) # add to tiles
        else:
            print(title)
            titles = np.append(titles,title)

        for line in lines:
            line = line.replace(' \n', '')

            h_line = line.replace(':', '').split(" ")
            cos_line = line.split(": ")

            res = [ h for h in hist_names if(h in h_line[0]) ]
            
            if bool(res) and "no_cuts" in h_line[0]:
                if title not in data_no_cuts:
                    data_no_cuts[title] = np.array([float(h_line[1]),float(h_line[2])])
                else:
                    data_no_cuts[title] = np.append(data_no_cuts[title],[float(h_line[1]),float(h_line[2])])
            elif bool(res) and "w_cuts" in h_line[0]:
                if title not in data_w_cuts:
                    data_w_cuts[title] = np.array([float(h_line[1]),float(h_line[2])])
                else:
                    data_w_cuts[title] = np.append(data_w_cuts[title],[float(h_line[1]),float(h_line[2])])
            elif "cos" in cos_line[0]:
                temp = np.array(cos_line[1].split(" ")).astype(float)

                if "%" in cos_line[0]:
                    if "no cut" in cos_line[0]:
                        if title not in cos_probs_no_cuts:
                            cos_probs_no_cuts[title] = temp
                        else:
                            cos_probs_no_cuts[title] = np.append(cos_probs_no_cuts[title],temp)
                    else:
                        if title not in cos_probs_w_cuts:
                            cos_probs_w_cuts[title] = temp
                        else:
                            cos_probs_w_cuts[title] = np.append(cos_probs_w_cuts[title],temp)
                elif "nums" in cos_line[0]:
                    temp_tot = [ temp[2*i] for i in range(0,int(len(temp)/2.0)) ]
                    temp_inv = [ temp[2*i+1] for i in range(0,int(len(temp)/2.0)) ]

                    if "no cut" in cos_line[0]:
                        if title not in cos_tots_no_cuts:
                            cos_tots_no_cuts[title] = temp_tot
                            cos_Se_Sa_no_cuts[title] = temp_inv
                        else:
                            cos_tots_no_cuts[title] = np.append(cos_tots_no_cuts[title],temp_tot)
                            cos_Se_Sa_no_cuts[title] = np.append(cos_Se_Sa_no_cuts[title],temp_inv)
                    else:
                        if title not in cos_tots_w_cuts:
                            cos_tots_w_cuts[title] = temp_tot
                            cos_Se_Sa_w_cuts[title] = temp_inv
                        else:
                            cos_tots_w_cuts[title] = np.append(cos_tots_w_cuts[title],temp_tot)
                            cos_Se_Sa_w_cuts[title] = np.append(cos_Se_Sa_w_cuts[title],temp_inv)
                elif "all" in cos_line[0]:
                    temp_tot = [ temp[2*i] for i in range(0,int(len(temp)/2.0)) ]
                    temp_inv = [ temp[2*i+1] for i in range(0,int(len(temp)/2.0)) ]

                    if title not in cos_all_s:
                        cos_all_s[title] = np.array([temp_tot,temp_inv])
                    else:
                        cos_all_s[title] = np.append(cos_all_s[title],np.array([temp_tot,temp_inv]))
                elif "minus" in cos_line[0]:
                    temp_tot = [ temp[2*i] for i in range(0,int(len(temp)/2.0)) ]
                    temp_inv = [ temp[2*i+1] for i in range(0,int(len(temp)/2.0)) ]

                    if title not in cos_h_minus:
                        cos_h_minus[title] = np.array([temp_tot,temp_inv])
                    else:
                        cos_h_minus[title] = np.append(cos_h_minus[title],np.array([temp_tot,temp_inv]))
                else:
                    cos_vals = temp

        in_file.close()

# ---------- Analysis ----------
print(titles)
# ----- Sort Input parameters
n_all = len(titles)
cos_vals_w_cuts = [ cos_val for cos_val in cos_vals if cos_val >= round(cos_min,1) and cos_val <= round(cos_max,1) ]

if type == "E":
    # E1 values
    E1_max = max(titles)
    E1_min = min(titles)
    E1_step = titles[1]-titles[0] # step size
    E1_plot = np.arange(E1_min,E1_max+1/4*E1_step,step=1/4*E1_step)
    E1_list = np.arange(0,n_all)

    if not np.all(titles[:-1] <= titles[1:]):
        print("E1 values not in order")
        title = titles.sort()

    p_4 = np.where(titles == 4)[0] # position for Belle 2 energy, 4 GeV

elif type == "contour":    
    n_m = len(m_values)
    n_g = len(g_values)

    # mA' values
    m_max = max(m_values)
    m_min = min(m_values)

    if n_m > 1:
        m_step = m_values[1]-m_values[0] # step size
        m_plot = np.arange(m_min,m_max+1/4*m_step,step=1/4*m_step)

    # g1' values
    g_max = max(g_values)
    g_min = min(g_values)

    if n_g > 1:
        g_step = g_values[1]-g_values[0] # step size
        g_plot = np.arange(g_min,g_max+1/4*g_step,step=1/4*g_step)

    # ----- Sort Input parameters
    m_list = np.arange(0,n_m)
    g_list = np.arange(0,n_g)

    if not np.all(m_values[:-1] <= m_values[1:]):
        print("m values not in order")

        m_values, m_list = (np.array(x) for x in zip(*sorted(zip(m_values,m_list), key=lambda pair:pair[0])))
        #titles_temp = titles.copy()
        # fix order of titles
        #for i_m in range(0,n_m):
            #titles[1 + i_m*n_g:1 + (i_m+1)*n_g] = [x for x in titles_temp[1 + m_list[i_m]*n_g:1 + (m_list[i_m]+1)*n_g] ]

    if not np.all(g_values[:-1] <= g_values[1:]):
        print("g values not in order")
    
        # fix order of titles
        #for i_m in range(0,n_m):
            #titles[1 + i_m*n_g:1 + (i_m+1)*n_g] = [x for _, x in sorted(zip(g_values,titles[1 + i_m*n_g:1 + (i_m+1)*n_g]), key=lambda pair: pair[0])]

        g_values, g_list = (np.array(x) for x in zip(*sorted(zip(g_values,g_list), key=lambda pair:pair[0])))

    # fix order of titles
    titles_temp = titles.copy()

    for title in titles_temp:
        a_title=title.split("_")
        ix1 = np.where(m_values == float(a_title[0]))[0][0]
        ix2 = np.where(g_values == float(a_title[1]))[0][0]
        titles[n_g*ix1+ix2]=title

    # ----- Move SM to correct position
    if n_m > 1:
        if ( m_min < 0 and m_max < 0 ): # both under zero
            m_plot_wZ = np.append(m_plot,[0.0]) # add zero to end
            m_wZ = np.append(m_values,[0.0])
        elif ( m_min > 0 and m_max > 0 ): # both over zero
            m_wZ = np.append([0.0],m_values)
            m_plot_wZ = np.append([0.0],m_plot) # add zero to beginning
        elif ( m_min < 0 and m_max > 0 ): # one above and one below
            m_plot_wZ = np.insert(m_plot,int(np.ceil(abs(m_min)/(1/4*m_step))),0.0)
            m_wZ = np.insert(m_values,int(np.ceil(abs(m_min)/(1/4*m_step))),0.0)
    if n_g > 1:
        if ( g_min < 0 and g_max < 0 ): # both under zero
            g_plot_wZ = np.append(g_plot,[0.0]) # add zero to end
            g_wZ = np.append(g_values,[0.0])
        elif ( g_min > 0 and g_max > 0 ): # both over zero
            g_plot_wZ = np.append([0.0],g_plot) # add zero to beginning
            g_wZ = np.append([0.0],g_values)
        elif ( g_min < 0 and g_max > 0 ): # one above and one below
            g_plot_wZ = np.insert(g_plot,int(np.ceil(abs(g_min)/(1/4*g_step))),0.0)
            g_wZ = np.insert(g_values,int(np.ceil(abs(g_min)/(1/4*g_step))),0.0)

    m_vals = []
    g_vals = []

    for j in range(0,n_g*n_m):
        m_vals.append(m_values[int(j/n_g)])
        g_vals.append(float(g_values[j%n_g]))    

    print(m_values)
    print(g_values)
"""
if type == "SM":
    # --- sum over data_no_cuts
    name_sum = "SUM"
    titles_wSUM = np.append(titles,name_sum)
    data_no_cuts[name_sum] = np.array([0,0]*n_h)

    for title in titles:
        if title != "All" and title != name_sum:
            data_no_cuts[name_sum] =  np.add(data_no_cuts[name_sum],data_no_cuts[title])
    n_all = len(titles)

    # -- sort titles
    titles = np.sort(titles)
    titles_wSUM = np.sort(titles_wSUM)
"""

# ----- Sort histogram values
h_no_cuts = {}
h_probs_no_cuts = {}
p_err_no_cuts = {}
cos_err_no_cuts = {}

h_w_cuts = {}
h_probs_w_cuts = {}
p_err_w_cuts = {}
cos_err_w_cuts = {}

# --- Data probabilites --- 
for i in range(0,n_h):
    h_no_cuts[hist_names[i]] = np.array([[data_no_cuts[title][2*i] for title in titles],[data_no_cuts[title][2*i+1] for title in titles]])
    h_w_cuts[hist_names[i]] = np.array([[data_w_cuts[title][2*i] for title in titles],[data_w_cuts[title][2*i+1] for title in titles]])

    # --- No cuts --- 
    # --- Probabilities
    tots = np.array([ (data_no_cuts[title][2*i] + data_no_cuts[title][2*i+1]) for title in titles ])
    list1 = [round(data_no_cuts[titles[j]][2*i]/tots[j],3) for j in range(0,n_all)]
    list2 = [round(data_no_cuts[titles[j]][2*i+1]/tots[j],3) for j in range(0,n_all)]
    h_probs_no_cuts[hist_names[i]] = np.array([list1,list2])

    # --- Errors
    list1 = [ percentage_error(data_no_cuts[titles[j]][2*i],tots[j]) for j in range(0,n_all)]
    list2 = [ percentage_error(data_no_cuts[titles[j]][2*i+1],tots[j]) for j in range(0,n_all)]
    p_err_no_cuts[hist_names[i]] = np.array([list1,list2])

    # --- With cuts --- 
    # --- Probabilities
    tots = np.array([ (data_w_cuts[title][2*i] + data_w_cuts[title][2*i+1]) for title in titles ])
    list1 = [round(data_w_cuts[titles[j]][2*i]/tots[j],3) for j in range(0,n_all)]
    list2 = [round(data_w_cuts[titles[j]][2*i+1]/tots[j],3) for j in range(0,n_all)]
    h_probs_w_cuts[hist_names[i]] = np.array([list1,list2])

    # --- Errors
    list1 = [ percentage_error(data_w_cuts[titles[j]][2*i],tots[j]) for j in range(0,n_all)]
    list2 = [ percentage_error(data_w_cuts[titles[j]][2*i+1],tots[j]) for j in range(0,n_all)]
    p_err_w_cuts[hist_names[i]] = np.array([list1,list2])

# ----- Crossing with 50% ----- 
cos_50 = {}

if type == "contour":
    cross_50 = [ np.where(np.round(cos_probs_w_cuts[title],1)==0.5)[0] for title in titles if g_const in title ]

    cos_SM_back = {}
    cos_SM_for = {}
    cross_SM_back = [ np.where(np.round(cos_probs_w_cuts[title],1)==round(SM_back,1))[0] for title in titles if g_const in title ]
    cross_SM_for = [ np.where(np.round(cos_probs_w_cuts[title],1)==round(SM_for,1))[0] for title in titles if g_const in title ]
else:
    cross_50 = [ np.where(np.round(cos_probs_no_cuts[title],1)==0.5)[0] for title in titles ]

iter = 0
for i in range(0,n_all):
    if type == "contour":
        if g_const in titles[i]:
            a_title = titles[i].split("_")
            #print(titles[i]," ",i)

            cos_50[a_title[0]] = np.array([ cos_vals[val] for val in cross_50[iter] ])
            cos_SM_back[a_title[0]] = np.array([ cos_vals[val] for val in cross_SM_back[iter] ])
            cos_SM_for[a_title[0]] = np.array([ cos_vals[val] for val in cross_SM_for[iter] ])

            iter += 1
    else:
        cos_50[titles[i]] = np.array([ cos_vals[val] for val in cross_50[iter] ])
        iter += 1

# ----- Cos value erros ----- 
if cos_tots_w_cuts:
    for title in titles:
        cos_err_w_cuts[title] = [ percentage_error(cos_Se_Sa_w_cuts[title][c_val],cos_tots_w_cuts[title][c_val]) for c_val in range(0,len(cos_vals)) ]
        cos_err_no_cuts[title] = [ percentage_error(cos_Se_Sa_no_cuts[title][c_val],cos_tots_no_cuts[title][c_val]) for c_val in range(0,len(cos_vals)) ]

# --- Structure for cos(theta) dependece ---
"""
test_n = len(cos_probs_w_cuts[titles[3]])
cos_a = np.sign( np.array([ cos_probs_w_cuts[titles[3]][i] - cos_probs_w_cuts[titles[3]][i+1] for i in range(0,test_n-1) if cos_probs_w_cuts[titles[3]][i+1] != 0.0 ])) # sign of accelerations
n_turns = np.array([ 1 for i in range(0,len(cos_a)-1) if cos_a[i] != cos_a[i+1] ]).sum()

if n_turns == 0.0:
    cross = round(cos_50[titles[3].split("_")[0]][0],1)
    cos_step = 0.1
    cross_m = cross
    cross_p = cross

    if cross < 0:
        cross_m = round(cos_min,1)
        cross_p = round(cos_max + cross,1)
    elif cross > 0:
        cross_m = round(cos_min + cross,1)
        cross_p = round(cos_max,1)
    print(cross," ",cross_p," ",cross_m)
    

    while cross_m >= cos_min and cross_p <= cos_max:
        cross_p = round(cross_p + 0.1,1)
        cross_m = round(cross_m - 0.1,1)
        #print(cross," ",cross_p," ",cross_m)

        temp1 = np.array([ cos_Se_Sa_w_cuts[titles[3]][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] >= round(cross_m,1) and cos_vals[c_val] <= round(cross_p,1) ])
        temp2 = np.array([ cos_tots_w_cuts[titles[3]][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] >= round(cross_m,1) and cos_vals[c_val] <= round(cross_p,1) ])
        per = float(temp1.sum())/float(temp2.sum())
        #print(per)
"""
#data_wo_zero = [x for x in cos_probs_w_cuts[titles[3]] if round(x,2) != 0]

#print(curve_fit(sigmod_func,cos_vals,cos_probs_w_cuts[titles[3]]))
"""
print(np.round(cos_probs_w_cuts[titles[3]],2))
print("---")
print(np.round([ np.cos(i) for i in np.linspace(0,np.pi/2,test_n) ],2))
print(np.round([ np.cos(i) for i in np.linspace(0,np.pi/2,test_n) ] - cos_probs_w_cuts[titles[3]],2))
print("---")
print(np.round([ 1/(1+np.exp(-i)) for i in np.linspace(2*np.pi,-2*np.pi,test_n) ],2))
print(np.round([ 1/(1+np.exp(-i)) for i in np.linspace(2*np.pi,-2*np.pi,test_n) ] - cos_probs_w_cuts[titles[3]],2))
print("---")
print(np.round([ math.erf(i) for i in np.linspace(np.pi/2,0,test_n) ],2))
print(np.round([ math.erf(i) for i in np.linspace(np.pi/2,0,test_n) ] - cos_probs_w_cuts[titles[3]],2))
"""
#elif n_turns == 1.0:

# --- Forwards/Backwards ---
for_no_cuts = np.array([])
back_no_cuts = np.array([])
for_w_cuts = np.array([])
back_w_cuts = np.array([])
for_no_cuts_errs = np.array([])
back_no_cuts_errs = np.array([])
for_w_cuts_errs = np.array([])
back_w_cuts_errs = np.array([])

for title in titles:
    # -- no cuts
    for_no_cuts,back_no_cuts,for_no_cuts_errs,back_no_cuts_errs=for_back_calcs(for_no_cuts,back_no_cuts,for_no_cuts_errs,back_no_cuts_errs,cos_Se_Sa_no_cuts[title],cos_tots_no_cuts[title],cos_vals,cos_min,cos_max)

    # -- w cuts
    for_w_cuts,back_w_cuts,for_w_cuts_errs,back_w_cuts_errs=for_back_calcs(for_w_cuts,back_w_cuts,for_w_cuts_errs,back_w_cuts_errs,cos_Se_Sa_w_cuts[title],cos_tots_w_cuts[title],cos_vals,cos_min,cos_max)
    """
    temp1 = np.array([ cos_Se_Sa_w_cuts[title][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] >= 0.0 and cos_vals[c_val] <= round(cos_max,1) ])
    temp2 = np.array([ cos_tots_w_cuts[title][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] >= 0.0 and cos_vals[c_val] <= round(cos_max,1) ])
    for_w_cuts = append_func_zero(for_w_cuts,temp1,temp2)

    # --- Errors
    for_errs = np.append(for_errs,percentage_error(temp1.sum(),temp2.sum()))

    # --- Backwards ---
    temp1 = np.array([ cos_Se_Sa_no_cuts[title][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] >= round(cos_min,1) and cos_vals[c_val] <= 0.0 ])
    temp2 = np.array([ cos_tots_no_cuts[title][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] >= round(cos_min,1) and cos_vals[c_val] <= 0.0 ])
    back_no_cuts = append_func_zero(back_no_cuts,temp1,temp2)

    temp1 = np.array([ cos_Se_Sa_w_cuts[title][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] >= round(cos_min,1) and cos_vals[c_val] <= 0.0 ])
    temp2 = np.array([ cos_tots_w_cuts[title][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] >= round(cos_min,1) and cos_vals[c_val] <= 0.0 ])
    back_w_cuts = append_func_zero(back_w_cuts,temp1,temp2)

    # --- Errors
    back_errs = np.append(back_errs,percentage_error(temp1.sum(),temp2.sum()))
    """

"""
if type == "E":
    # ---------- Fitting ----------
    fit_funcs = {}
    
    for i in range(1,n_h):
        fit0, cov0 = curve_fit(f1,titles,h_probs_no_cuts[hist_names[i]][0])
        fit1, cov1 = curve_fit(f1,titles,h_probs_no_cuts[hist_names[i]][1])
        fit_funcs[hist_names[i]] = np.array([fit0,fit1])
"""

# --- Arrange data ---
if type == "contour":
    # structure for contour plot
    z_g = [] # as a function of mass
    z_g_err = []
    z_m = [] # as a function of coupling constant
    z_m_err = []

    for_g = []
    back_g = []
    for_m = []
    back_m = []

    for i in range(0,n_all):
        if i < n_m: # constant mass, vary coupling
            z_m.append([ h_probs_no_cuts[hist_names[2]][1][g + (i-1)*n_g] for g in range(0,n_g) ])
            z_m_err.append([ p_err_no_cuts[hist_names[2]][1][g + (i-1)*n_g] for g in range(0,n_g) ])

            for_m.append([ for_no_cuts[g + (i-1)*n_g] for g in range(0,n_g) ])
            back_m.append([ back_no_cuts[g + (i-1)*n_g] for g in range(0,n_g) ])
        if i < n_g: # constant coupling, vary mass
            z_g.append([ h_probs_no_cuts[hist_names[2]][1][i + m*n_g] for m in range(0,n_m) ])
            z_g_err.append([ p_err_no_cuts[hist_names[2]][1][i + m*n_g] for m in range(0,n_m) ])

            for_g.append([ for_no_cuts[i + m*n_g] for m in range(0,n_m) ])
            back_g.append([ back_no_cuts[i + m*n_g] for m in range(0,n_m) ])

# ---------- Printing ----------
"""
if type == "SM":
    print('all')
    print("w cuts")
    print("sum: ", round(float(np.array(cos_Se_Sa_w_cuts['all']).sum())/float(np.array(cos_tots_w_cuts['all']).sum()),3))
    print("average: ", round(float(cos_probs_w_cuts['all'].sum())/float(len(cos_probs_w_cuts['all'])),3))

    print("no cuts")
    print("sum: ", round(float(np.array(cos_Se_Sa_no_cuts['all']).sum())/float(np.array(cos_tots_no_cuts['all']).sum()),3))
    print("average: ", round(float(cos_probs_no_cuts['all'].sum())/float(len(cos_probs_w_cuts['all'])),3))
    temp1 = np.array([ cos_Se_Sa_no_cuts['all'][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] > round(cos_min,1) and cos_vals[c_val] < round(cos_max,1) ])
    temp2 = np.array([ cos_tots_no_cuts['all'][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] > round(cos_min,1) and cos_vals[c_val] < round(cos_max,1) ])
    print(round(float(temp1.sum())/float(temp2.sum()),3))

    print("from analysis")
    print("no cuts: ", h_probs_no_cuts['hsE1'][1][0])
    print("no cuts err: ", p_err_no_cuts['hsE1'][1][0])
    print("w cuts: ", h_probs_w_cuts['hsE1'][1][0])
    print("w cuts err: ", p_err_w_cuts['hsE1'][1][0])

    print("forwards/backwards")
    print("w. cuts")
    print("forward: ",round(for_w_cuts[0],3))
    print("forward err: ",round(for_w_cuts_errs[0],3))
    print("backward: ",round(back_w_cuts[0],3))
    print("backward err: ",round(back_w_cuts_errs[0],3))
"""

# ---------- Plotting ----------
# --- plot 1: Percentages
fig1a, ax1a = plt.subplots(1,2)
fig1b, ax1b = plt.subplots(1,2)

for i in range(2,n_h):
    hist_name = hist_names[i].replace("hs","").replace("E1","$S(E)=1$").replace("Em1","$S(E)=-1$")

    for ax in [ax1a,ax1b]:
        if ax in np.array([ax1a]):
            prob = h_probs_no_cuts
            h = h_no_cuts
        elif ax in [ax1b]:
            prob = h_probs_w_cuts
            h = h_w_cuts

        ax[i-2].scatter(titles,prob[hist_names[i]][0],label="$S(A)=-1$",marker=".")
        ax[i-2].scatter(titles,prob[hist_names[i]][1],label="$S(A)=1$",marker=".")

        if type == "E":
            """
            fit0=fit_funcs[hist_names[i]][0]
            y0 = f1(E1_plot,fit0[0],fit0[1],fit0[2],fit0[3],fit0[4],fit0[5],fit0[6])
            ax1a[i-2].plot(E1_plot, y0, linewidth=0.8, alpha=0.8)

            fit1=fit_funcs[hist_names[i]][1]
            y1 = f1(E1_plot,fit1[0],fit1[1],fit1[2],fit1[3],fit1[4],fit1[5],fit1[6])
            ax1a[i-2].plot(E1_plot,y1, linewidth=0.8, alpha=0.8)

            # Intersection
            idx = np.argwhere(np.diff(np.sign(y0 - y1))).flatten()
            ax1a[i-2].scatter(E1_plot[idx],y0[idx],color="grey")
            ax1a[i-2].vlines(E1_plot[idx],0.0,1.0,color="grey",ls="-",label="Intersection: "+str(float(E1_plot[idx]))+" GeV", linewidth=0.8, alpha=0.8)
            """
            # Belle 2 values
            ax1a[i-2].scatter([4],prob[hist_names[i]][0][p_4],color='r',marker=".")
            ax1a[i-2].scatter([4],prob[hist_names[i]][1][p_4],color='r',marker=".")
            ax1a[i-2].vlines(4,0.0,1.0,color="red",ls="-",label="Belle 2: 4 GeV", linewidth=0.8, alpha=0.8)

            ax1a[i-2].set_xlabel("$E_1$ [GeV]")
            ax1a[i-2].set_xticks(np.round(np.linspace(0,E1_max,6),0))

            # remove other axes
            ax1a[i-2].spines['left'].set_position('zero')

        elif type == "SM":
            ax[i-2].tick_params(axis='x', rotation=90)

            l_ones = np.ones(len(prob[hist_names[i]][0]))
        
            Nerr0 = sqrt(np.divide(l_ones,h[hist_names[i]][0])*(l_ones + np.divide(l_ones,prob[hist_names[i]][0])))
            Nerr1 = sqrt(np.divide(l_ones,h[hist_names[i]][1])*(l_ones + np.divide(l_ones,prob[hist_names[i]][1])))

            ax[i-2].errorbar(titles,prob[hist_names[i]][0],yerr=prob[hist_names[i]][0]*Nerr0,ls="")
            ax[i-2].errorbar(titles,prob[hist_names[i]][1],yerr=prob[hist_names[i]][1]*Nerr1,ls="")

        ax[i-2].spines['right'].set_color('none')
        ax[i-2].spines['top'].set_color('none')
        #ax[i-2].spines['bottom'].set_position('zero')

        ax[i-2].set_yticks(np.linspace(0,1.0,11))

        ax[i-2].set_title(hist_name)
        ax[i-2].grid(True)

        ax[1].legend()
        ax[0].set_ylabel(prob_text)

fig1a.suptitle(prob_text+" of $s(\\gamma)=s(e^-)$")
fig1b.suptitle(prob_text+" of $s(\\gamma)=s(e^-)$")
fig1a.savefig("probs_no_cuts.jpg" , bbox_inches='tight', dpi=250)
fig1b.savefig("probs_w_cuts.jpg" , bbox_inches='tight', dpi=250)


if type == "contour":
    [X, Y] = np.meshgrid(m_values,g_values)

    # --- plot 2: Contour Levels
    fig2, ax2 = plt.subplots(1,1)
    z_step = (np.max(z_g)-np.min(z_g))/100.

    if n_g > 1 and n_m > 1:
        levels = np.arange(np.min(z_g), np.max(z_g)+z_step, z_step)
        im = ax2.contourf(X,Y,z_g, levels=levels)
        fig2.colorbar(im) #, label="$S(A)=1 & S(E)=1 "+prob_text)

    # --- plot 3 & 4: SE1 and SA1 as a function of mass and coupling
    fig3, ax3 = plt.subplots(1,1)
    fig4, ax4 = plt.subplots(1,1)

    iter=0
    for i in range(0,n_all):
        if i < n_g: # loop over coupling constants
            ax3.scatter(m_values,z_g[i], label=str("{:.1e}".format(g_values[i-1])),marker=".")
            ax3.errorbar(m_values,z_g[i],yerr=z_g_err[i],ls="")

        if i < n_m: # loop over masses
            ax4.scatter(g_values,z_m[i], label=str(m_values[i-1]),marker=".")
            ax4.errorbar(g_values,z_m[i],yerr=z_m_err[i],ls="")

    # ----- set title, axses, etc. 
    fig_axes = [ax2,ax3,ax4]
    figs = [fig2,fig3,fig4]

    plot_titles = ["contour","SE1_SA1_m","SE1_SA1_g"]

    for ax in fig_axes:
        ax.set_title("S(A)=1 & S(E)=1 "+prob_text)

        # --- set x-axis
        if ax in [ax2,ax3]:
            ax.set_xlabel(m_text+" [GeV]")
            
        elif ax in [ax4]:
            ax.set_xlabel(g_text)
            ax.set_xscale("log")
        
        # --- set y-axis
        if ax in [ax2]:
            ax.set_ylabel(g_text)

            if proc == "ALP":
                lab = ax.get_yticks()
                ax.set_yticklabels([np.format_float_scientific(lab_el, unique=False, precision=1) for lab_el in lab])
                ax.set_yscale("log")

        elif ax in [ax3,ax4]:
            ax.set_ylabel(prob_text)
        
        # --- set legends
        if ax in [ax4]:
            figs[iter].legend(bbox_to_anchor=(0.9,1), loc="upper left", ncol=1,title=m_text+" [GeV]")
        elif ax in [ax3]:
            figs[iter].legend(bbox_to_anchor=(0.9,1), loc="upper left", ncol=1,title=g_text)

        # --- set grid
        ax.grid()

        # --- save plot
        figs[iter].savefig(plot_titles[iter]+".jpg" , bbox_inches='tight', dpi=250)

        iter+=1

if bool(cos_probs_w_cuts): # if there is data_no_cuts in cos
    if type == "contour": 
        # --- plot 5: angular distribution data_no_cuts for mass and coupling
        fig5, ax5 = plt.subplots(1,1)#2) # no cuts
        fig6, ax6 = plt.subplots(1,1)#2) # no cuts
        fig7, ax7 = plt.subplots(1,1)#2) # w cuts
        fig8, ax8 = plt.subplots(1,1)#2) # w cuts

        for title in titles:
            a_title=np.array(title.split("_"))
            if g_const in a_title[1] and ".0" in str(a_title[0]): # constant coupling (1*10^(-3)), vary mass
                #ax5[0].scatter(cos_vals,cos_probs[title], label=title)
                #ax5[1].scatter(cos_vals,cos_probs_no_cuts[title], label=title.replace("_"+g_const,""))
                label_text = str(title.replace("_"+g_const,""))
                ax5.scatter(cos_vals,cos_probs_no_cuts[title], label=label_text,marker=".")
                ax7.scatter(cos_vals,cos_probs_w_cuts[title], label=label_text,marker=".")

                """
                if proc == "DP":
                    print(title)
                    fit = curve_fit(pol_cos,cos_vals,cos_probs_no_cuts[title])[0]
                    cos_plot = np.linspace(min(cos_vals),max(cos_vals),100)
                    ax5.plot(cos_plot,pol_cos(cos_plot,*fit))
                """

            if m_const in a_title[0]: # constant coupling (1 GeV), vary mass
                #ax6[0].scatter(cos_vals,np.round(cos_probs[title],2), label=title)
                #ax6[1].scatter(cos_vals,np.round(cos_probs_no_cuts[title],2), label=str("{:.1e}".format(float(title.replace(m_const+"_","")))))
                label_text = str("{:.1e}".format(float(title.replace(m_const+"_",""))))
                ax6.scatter(cos_vals,np.round(cos_probs_no_cuts[title],2), label=label_text,marker=".")
                ax8.scatter(cos_vals,np.round(cos_probs_w_cuts[title],2), label=label_text,marker=".")

                """
                if proc == "DP":
                    fit = curve_fit(pol_cos,cos_vals,cos_probs_no_cuts[title])[0]
                    cos_plot = np.linspace(min(cos_vals),max(cos_vals),100)
                    ax6.plot(cos_plot,pol_cos(cos_plot,*fit))
                """

        for ax in [ax5,ax6,ax7,ax8]: #[ax5[0],ax5[1],ax6[0],ax6[1]]:
            ax.set_ylim([0.0,1.0])
            #ax.set_yticks(np.linspace(0.0,1.0,11))

            ax.grid(True)
            ax.grid(True)

            if ax in [ax5,ax7]:
                ax.legend(title=m_text+" [GeV]")
            elif ax in [ax6,ax8]:
                ax.legend(title=g_text)
            ax.set_ylabel(prob_text)
            ax.set_xlabel("$cos(\\theta_{\\gamma})$")
            ax.axvline(0,color='black', linewidth=0.5, alpha=0.8) # add line at x = 0

        fig5.suptitle("$s(\\gamma) = s(e^-)$ as a function of $cos(\\theta_{\\gamma})$ and "+m_text)
        ax5.set_title(g_text+"= "+str("{:.1e}".format(float(g_const))))
        fig7.suptitle("$s(\\gamma) = s(e^-)$ as a function of $cos(\\theta_{\\gamma})$ and "+m_text)
        ax7.set_title(g_text+"= "+str("{:.1e}".format(float(g_const))))

        fig6.suptitle("$s(\\gamma) = s(e^-)$ as a function of $cos(\\theta_{\\gamma})$ and "+g_text)
        ax6.set_title(m_text+"= "+m_const+" GeV")
        fig8.suptitle("$s(\\gamma) = s(e^-)$ as a function of $cos(\\theta_{\\gamma})$ and "+g_text)
        ax8.set_title(m_text+"= "+m_const+" GeV")

        fig5.savefig("cos_m_no_cuts.jpg" , bbox_inches='tight', dpi=250)
        fig7.savefig("cos_m_w_cuts.jpg" , bbox_inches='tight', dpi=250)
        fig6.savefig("cos_g_no_cuts.jpg" , bbox_inches='tight', dpi=250)
        fig8.savefig("cos_g_w_cuts.jpg" , bbox_inches='tight', dpi=250)

    elif type == "SM" and SM_more == 1:
        # --- plot 5: angular distribution data_no_cuts
        fig5, ax5 = plt.subplots(2,2)
        fig6, ax6 = plt.subplots(2,2)

        for title in titles:
            title_text=title
            if "e" in title:
                ix=0
                iy=0
            elif "nu" in title:
                ix=0
                iy=1
            elif "all" in title:
                ix=1
                iy=0
                title_text="Total SM Bg"
            else:
                ix=1
                iy=1

            i_cut = np.where(cos_probs_no_cuts[title]==0)[0]
            x_vals= [cos_vals[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut]

            y_vals = [ cos_probs_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
            ax5[ix][iy].scatter(x_vals,y_vals, label=title_text, marker=".")

            y_vals = [ cos_probs_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
            ax6[ix][iy].scatter(x_vals,y_vals, label=title_text, marker=".")

            if cos_probs_w_cuts:
                y_vals = [ cos_probs_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                y_errs = [ cos_err_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                ax5[ix][iy].errorbar(x_vals,y_vals,yerr=y_errs,ls="")

                y_vals = [ cos_probs_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                y_errs = [ cos_err_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                ax6[ix][iy].errorbar(x_vals,y_vals,yerr=y_errs,ls="")

            """ 
            ax5[ix][iy].scatter(cos_vals,cos_probs_no_cuts[title], label=title_text, marker=".")
            ax6[ix][iy].scatter(cos_vals,cos_probs_w_cuts[title], label=title_text, marker=".")
            if cos_probs_w_cuts:
                ax5[ix][iy].errorbar(cos_vals,cos_probs_no_cuts[title],yerr=cos_err_no_cuts[title],ls="")
                ax6[ix][iy].errorbar(cos_vals,cos_probs_w_cuts[title],yerr=cos_err_w_cuts[title],ls="")
            """

        for ax in [ax5,ax6]:
            for i in [0,1]:
                ax[i][0].set_ylabel(prob_text)
                ax[i][1].set_yticklabels([])
                ax[i][0].set_yticks(np.linspace(0.0,1.0,11),["0","10","20","30","40","50","60","70","80","90","100"])


                ax[1][i].set_xlabel("$cos(\\theta_{\\gamma})$")
                ax[1][i].set_xlim([-1.1,1.1])
                ax[0][i].set_xticklabels([])

                for j in [0,1]:
                    ax[i][j].legend()
                    ax[i][j].grid(True)
                    ax[i][j].set_yticks(np.linspace(0,1.0,11))
                    ax[i][j].axvline(0,color='black', linewidth=0.5, alpha=0.8) # add line at x = 0

        fig5.subplots_adjust(wspace=0.05,hspace=0.1)
        fig6.subplots_adjust(wspace=0.05,hspace=0.1)
        #fig5.suptitle("$s(\\gamma) = s(e^-)$ for angular distribution of $\\gamma$")
        #fig6.suptitle("$s(\\gamma) = s(e^-)$ for angular distribution of $\\gamma$")
        fig5.savefig("cos_no_cuts.jpg" , bbox_inches='tight', dpi=250)
        fig6.savefig("cos_w_cuts.jpg" , bbox_inches='tight', dpi=250)

    else:
        # --- plot 5: angular distribution data_no_cuts
        fig5, ax5 = plt.subplots(1,1)#2)
        fig6, ax6 = plt.subplots(1,1)#2)

        for title in titles:
            i_cut = np.where(cos_probs_no_cuts[title]==0)[0]
            x_vals= [cos_vals[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut]

            title_text=str(title)
            #ax5.scatter(cos_vals,cos_probs_no_cuts[title], label=title_text,marker=".")
            #ax6.scatter(cos_vals,cos_probs_w_cuts[title], label=title_text,marker=".")

            y_vals = [ cos_probs_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
            ax5.scatter(x_vals,y_vals, label=title_text,marker=".")

            y_vals = [ cos_probs_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
            ax6.scatter(x_vals,y_vals, label=title_text,marker=".")

            if cos_probs_w_cuts:
                #ax5.errorbar(cos_vals,cos_probs_no_cuts[title],yerr=cos_err_no_cuts[title],ls="")
                #ax6.errorbar(cos_vals,cos_probs_w_cuts[title],yerr=cos_err_w_cuts[title],ls="")

                y_vals = [ cos_probs_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                y_errs = [ cos_err_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                ax5.errorbar(x_vals,y_vals,yerr=y_errs,ls="")

                y_vals = [ cos_probs_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                y_errs = [ cos_err_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                ax6.errorbar(x_vals,y_vals,yerr=y_errs,ls="")

        for ax in [ax5,ax6]:
            ax.set_ylabel(prob_text)
            ax.set_yticks(np.linspace(0,1.0,11))
            ax.legend()
            ax.grid(True)
            ax.set_xlabel("$cos(\\theta_{\\gamma})$")
            ax.set_xlim(-1,1)
            ax.axvline(0,color='black', linewidth=0.5, alpha=0.8) # add line at x = 0

        if type == "E":
            fig5.suptitle("Variation of positron beam energy")
            fig6.suptitle("Variation of positron beam energy")
        elif type == "SM":
            fig5.suptitle(prob_text+" of $s(\\gamma)=s(e^-)$ func of $cos(\\theta_{\\gamma})$")
            fig6.suptitle(prob_text+" of $s(\\gamma)=s(e^-)$ func of $cos(\\theta_{\\gamma})$")

        fig5.savefig("cos_no_cuts.jpg" , bbox_inches='tight', dpi=250)
        fig6.savefig("cos_w_cuts.jpg" , bbox_inches='tight', dpi=250)

# --- plot 9: Forwards/Backwards
if type == "SM":
    fig9, ax9 = plt.subplots(1,2)
    ax9[0].scatter(titles,for_w_cuts,marker=".")
    ax9[0].errorbar(titles,for_w_cuts,yerr=for_w_cuts_errs,ls="")

    ax9[1].scatter(titles,back_w_cuts,marker=".")
    ax9[1].errorbar(titles,back_w_cuts,yerr=back_w_cuts_errs,ls="")

    for i in [0,1]:
        ax9[i].set_yticks(np.linspace(0,1.0,11))
        ax9[i].tick_params(axis='x', rotation=90)
        ax9[i].grid(True)
    ax9[0].set_ylabel(prob_text)

    ax9[0].set_title("$cos(\\theta_\\gamma) > 0$")
    ax9[1].set_title("$cos(\\theta_\\gamma) < 0$")

    fig9.suptitle("$s(\\gamma) = s(e^-)$ "+prob_text)
    fig9.savefig("for_back.jpg", bbox_inches='tight', dpi=250)

elif type == "contour":
    # --- plot 9: Contour Levels
    fig9a, ax9a = plt.subplots(1,1)
    fig9b, ax9b = plt.subplots(1,1)
    
    [X, Y] = np.meshgrid(m_values,g_values)
    for_step = (np.max(for_g)-np.min(for_g))/100.
    back_step = (np.max(back_g)-np.min(back_g))/100.

    # --- plot 1: function of mass
    fig10a, ax10a = plt.subplots(1,1)
    fig10b, ax10b = plt.subplots(1,1)

    if n_g > 1 and n_m > 1:
        levels = np.arange(np.min(for_g), np.max(for_g)+for_step, for_step)
        ima = ax9a.contourf(X,Y,for_g, levels=levels)
        fig9a.colorbar(ima)

        levels = np.arange(np.min(back_g), np.max(back_g)+back_step, back_step)
        imb = ax9b.contourf(X,Y,back_g, levels=levels)
        fig9b.colorbar(imb)
        
        for ax in [ax9a,ax9b,ax10a,ax10b]:
            # --- set x-axis
            ax.set_xlabel(m_text+" [GeV]")

            # --- set grid
            ax.grid()
        
            # --- set y-axis
            if ax in [ax9a,ax9b]:
                ax.set_ylabel(g_text)

                if proc == "ALP":
                    lab = ax.get_yticks()
                    ax.set_yticklabels([np.format_float_scientific(lab_el, unique=False, precision=1) for lab_el in lab])
                    ax.set_yscale("log")
            elif ax in [ax10a,ax10b]:
                ax.set_ylabel(prob_text)

            if ax in [ax9a,ax10a]:
                ax.set_title("$cos(\\theta_\\gamma) > 0$")
            elif ax in [ax9b,ax10b]:
                ax.set_title("$cos(\\theta_\\gamma) < 0$")

    for fig in [fig9a,fig9b,fig10a,fig10b]:
        fig.suptitle("$s(\\gamma) = s(e^-)$ "+prob_text)

    for i in range(0,n_all):
        if i > n_g and i > n_m:
            break
        else:
            if n_m > 1 and i < n_g: # loop over coupling constants
                ax10a.scatter(m_values,for_g[i], label=str("{:.1e}".format(g_values[i])),marker=".")
                ax10b.scatter(m_values,back_g[i], label=str("{:.1e}".format(g_values[i])),marker=".")

    ax10a.plot(m_values,[SM_for]*n_m, label="SM value",color="r")
    ax10b.plot(m_values,[SM_back]*n_m, label="SM value",color="r")

    fig10a.legend(loc='upper left',bbox_to_anchor=(0.125, 0.85), ncol=2,title=g_text)
    fig10b.legend(loc='lower left',bbox_to_anchor=(0.125, 0.15), ncol=2,title=g_text)


    fig9a.savefig("for.jpg", bbox_inches='tight', dpi=250)
    fig9b.savefig("back.jpg", bbox_inches='tight', dpi=250)
    fig10a.savefig("for_m.jpg", bbox_inches='tight', dpi=250)
    fig10b.savefig("back_m.jpg", bbox_inches='tight', dpi=250)
