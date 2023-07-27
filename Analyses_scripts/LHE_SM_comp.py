# ---------- Libaries ----------
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('tableau-colorblind10')

import sys
import os
import math
from math import cos, acos, sin
import numpy as np
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.special import expit
import re

from sympy.parsing.mathematica import parse_mathematica
from sympy.parsing.mathematica import mathematica
import scipy.integrate as integrate
from scipy.integrate import trapz
from sympy import var
from sympy import N

from Extra.stat_extra import *
from Extra.fit_funcs import *
from Extra.DP_ALP_extra import *

# ---------- Variables ----------
direc = sys.argv[1] # directory path
prob_text = "Percentage (%)" #"Probability"

cl = 0.95   # 95% condidence level
dof = 1  # degrees of freedom per bin
chi_crit_1 = stats.chi2.ppf(cl, dof) # critical chi-squared value for 1 dof

# ---------- Mass and coupling variables
m_values = np.array([])
g_values = np.array([])

g_text = "$g$"
m_text = "$m$"

if ( len(sys.argv) < 3 ): # less than three command line arguments
    print("Not correct number of command line arguments")
    quit()
else:
    paper_on = int(sys.argv[3]) # plots for paper (1 for on, 2 for poster)

    proc = sys.argv[2] # type of process
    if proc == "DP":
        g_text = "$g_{\gamma'}$"
        g_const = "0.001"
        m_text = "$M_{\gamma'}$"
    elif "ALP" in proc:
        g_text = "$f_{a}$"
        g_const = "40000.0"
        m_text = "$M_{ALP}$"
    elif proc == "scalar":
        g_text = "g_i"
        g_const = ""
        m_text = "$M_X$"
        
    m_const = "1.0"

Eng1=4
Eng2=7

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
cos_min = np.cos(170*(np.pi/180))
cos_max = np.cos(10*(np.pi/180))

pattern = r'[0-9]'

# ---------- Read through files ----------
# return all files as a list
for file in os.listdir(direc):
    if file.endswith(ending) and "out" not in file:
        #print(file)
        in_file = open(os.path.join(direc, file),'r')
        lines = in_file.readlines()

        title = file.replace(ending,'').replace("E_theta_","")
    
        if "SM" in title:
            title = "SM"
        else:
            if paper_on == 2:
                if "ALPb" in title:
                    title_temp = "ALP t/u-channel"
                    if "3" in title:
                        title = title_temp + " 3 GeV"
                    elif "9" in title:
                        title = title_temp + " 9 GeV"
                    else:
                        title = title_temp
                elif "ALPs_DP" in title:
                    title_temp = "ALP s-channel & \n Dark Photon"
                    if "3" in title:
                        title = title_temp + " 3 GeV"
                    else:
                        title = title_temp
                elif "DP" in title and "9" in title:
                    title = "Dark Photon 9 GeV"

                if m_values.size == 0:
                    m_values = np.array([9])
                    g_values = np.array([9])
            else:
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

if "SM" not in titles:
    print("SM file not found")
    quit()
else:
    SM_pos = np.where(titles=="SM")

# ---------- Analysis ----------
# ----- Sort Input parameters
n_all = len(titles)
cos_vals_w_cuts = [ cos_val for cos_val in cos_vals if cos_val >= round(cos_min,1) and cos_val <= round(cos_max,1) ]

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

if not np.all(g_values[:-1] <= g_values[1:]):
    print("g values not in order")
    g_values, g_list = (np.array(x) for x in zip(*sorted(zip(g_values,g_list), key=lambda pair:pair[0])))

# fix order of titles
titles_temp = titles.copy()
titles_temp = np.delete(titles_temp,SM_pos)

if proc != "all":
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

#print(m_values)
#print(g_values)
print(titles)

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
    h_no_cuts[hist_names[i]] = np.array([[data_no_cuts[title][2*i] for title in titles if title!="SM"],[data_no_cuts[title][2*i+1] for title in titles if title!="SM"]])
    h_w_cuts[hist_names[i]] = np.array([[data_w_cuts[title][2*i] for title in titles if title!="SM"],[data_w_cuts[title][2*i+1] for title in titles if title!="SM"]])

    # --- No cuts --- 
    # --- Probabilities
    tots = np.array([ (data_no_cuts[title][2*i] + data_no_cuts[title][2*i+1]) for title in titles if title!="SM"])
    list1 = [round(data_no_cuts[titles[j]][2*i]/tots[j],3) for j in range(0,n_all-1) if j!= SM_pos]
    list2 = [round(data_no_cuts[titles[j]][2*i+1]/tots[j],3) for j in range(0,n_all-1) if j!= SM_pos]
    h_probs_no_cuts[hist_names[i]] = np.array([list1,list2])

    # --- Errors
    list1 = [ percentage_error(data_no_cuts[titles[j]][2*i],tots[j]) for j in range(0,n_all-1) if j!= SM_pos]
    list2 = [ percentage_error(data_no_cuts[titles[j]][2*i+1],tots[j]) for j in range(0,n_all-1) if j!= SM_pos]
    p_err_no_cuts[hist_names[i]] = np.array([list1,list2])

    # --- With cuts --- 
    # --- Probabilities
    tots = np.array([ (data_w_cuts[title][2*i] + data_w_cuts[title][2*i+1]) for title in titles if title!="SM"])
    list1 = [round(data_w_cuts[titles[j]][2*i]/tots[j],3) for j in range(0,n_all-1) if j!= SM_pos]
    list2 = [round(data_w_cuts[titles[j]][2*i+1]/tots[j],3) for j in range(0,n_all-1) if j!= SM_pos]
    h_probs_w_cuts[hist_names[i]] = np.array([list1,list2])

    # --- Errors
    list1 = [ percentage_error(data_w_cuts[titles[j]][2*i],tots[j]) for j in range(0,n_all-1) if j!= SM_pos]
    list2 = [ percentage_error(data_w_cuts[titles[j]][2*i+1],tots[j]) for j in range(0,n_all-1) if j!= SM_pos]
    p_err_w_cuts[hist_names[i]] = np.array([list1,list2])

# ----- Crossing with 50% ----- 
if proc != "all":
    cos_50 = {}
    cross_50 = [ np.where(np.round(cos_probs_w_cuts[title],1)==0.5)[0] for title in titles if g_const in title ]

    cos_SM_back = {}
    cos_SM_for = {}
    #cross_SM_back = [ np.where(np.round(cos_probs_w_cuts[title],1)==round(SM_back,1))[0] for title in titles if g_const in title ]
    #cross_SM_for = [ np.where(np.round(cos_probs_w_cuts[title],1)==round(SM_for,1))[0] for title in titles if g_const in title ]

    iter = 0
    for i in range(0,n_all):
        if g_const in titles[i]:
            a_title = titles[i].split("_")
            #print(titles[i]," ",i," ",iter)

            cos_50[a_title[0]] = np.array([ cos_vals[val] for val in cross_50[iter] ])
            #cos_SM_back[a_title[0]] = np.array([ cos_vals[val] for val in cross_SM_back[iter] ])
            #cos_SM_for[a_title[0]] = np.array([ cos_vals[val] for val in cross_SM_for[iter] ])
            iter += 1

# ----- Cos value erros ----- 
if cos_tots_w_cuts:
    for title in titles:
        cos_err_w_cuts[title] = [ percentage_error(cos_Se_Sa_w_cuts[title][c_val],cos_tots_w_cuts[title][c_val]) for c_val in range(0,len(cos_vals)) ]
        cos_err_no_cuts[title] = [ percentage_error(cos_Se_Sa_no_cuts[title][c_val],cos_tots_no_cuts[title][c_val]) for c_val in range(0,len(cos_vals)) ]
    cos_err_w_cuts["SM"] = [ percentage_error(cos_Se_Sa_w_cuts["SM"][c_val],cos_tots_w_cuts["SM"][c_val]) for c_val in range(0,len(cos_vals)) ]
    cos_err_no_cuts["SM"] = [ percentage_error(cos_Se_Sa_no_cuts["SM"][c_val],cos_tots_no_cuts["SM"][c_val]) for c_val in range(0,len(cos_vals)) ]

# --- Structure for cos(theta) dependece ---
if proc != "all":
    test_n = len(cos_probs_w_cuts[titles[3]])
    cos_a = np.sign( np.array([ cos_probs_w_cuts[titles[3]][i] - cos_probs_w_cuts[titles[3]][i+1] for i in range(0,test_n-1) if cos_probs_w_cuts[titles[3]][i+1] != 0.0 ])) # sign of accelerations
    n_turns = np.array([ 1 for i in range(0,len(cos_a)-1) if cos_a[i] != cos_a[i+1] ]).sum()

    if n_turns == 0.0:
        if len(cos_50[titles[3].split("_")[0]]) != 0:
            cross = round(cos_50[titles[3].split("_")[0]][0],1)
            cos_step = 0.1
            cross_m = cross
            cross_p = cross

            while cross_m >= cos_min and cross_p <= cos_max:
                cross_p = round(cross_p + 0.1,1)
                cross_m = round(cross_m - 0.1,1)

                temp1 = np.array([ cos_Se_Sa_w_cuts[titles[3]][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] >= round(cross_m,1) and cos_vals[c_val] <= round(cross_p,1) ])
                temp2 = np.array([ cos_tots_w_cuts[titles[3]][c_val] for c_val in range(0,len(cos_vals)) if cos_vals[c_val] >= round(cross_m,1) and cos_vals[c_val] <= round(cross_p,1) ])
                per = float(temp1.sum())/float(temp2.sum())

# --- Forwards/Backwards ---
for_no_cuts = np.array([])
back_no_cuts = np.array([])
for_w_cuts = np.array([])
back_w_cuts = np.array([])
for_no_cuts_errs = np.array([])
back_no_cuts_errs = np.array([])
for_w_cuts_errs = np.array([])
back_w_cuts_errs = np.array([])

# -- such that SM is zero position
for_no_cuts,back_no_cuts,for_no_cuts_errs,back_no_cuts_errs=for_back_calcs(for_no_cuts,back_no_cuts,for_no_cuts_errs,back_no_cuts_errs,cos_Se_Sa_no_cuts["SM"],cos_tots_no_cuts["SM"],cos_vals,cos_min,cos_max)
for_w_cuts,back_w_cuts,for_w_cuts_errs,back_w_cuts_errs=for_back_calcs(for_w_cuts,back_w_cuts,for_w_cuts_errs,back_w_cuts_errs,cos_Se_Sa_w_cuts["SM"],cos_tots_w_cuts["SM"],cos_vals,cos_min,cos_max)

for title in titles:
    if title != "SM":
        # -- no cuts
        for_no_cuts,back_no_cuts,for_no_cuts_errs,back_no_cuts_errs=for_back_calcs(for_no_cuts,back_no_cuts,for_no_cuts_errs,back_no_cuts_errs,cos_Se_Sa_no_cuts[title],cos_tots_no_cuts[title],cos_vals,cos_min,cos_max)
        # -- w cuts
        for_w_cuts,back_w_cuts,for_w_cuts_errs,back_w_cuts_errs=for_back_calcs(for_w_cuts,back_w_cuts,for_w_cuts_errs,back_w_cuts_errs,cos_Se_Sa_w_cuts[title],cos_tots_w_cuts[title],cos_vals,cos_min,cos_max)

# --- Arrange data ---
# structure for contour plot
z_g_no_cuts = [] # as a function of mass
z_g_err_no_cuts = []
z_m_no_cuts = [] # as a function of coupling constant
z_m_err_no_cuts = []

for_g_no_cuts = []
back_g_no_cuts = []
for_m_no_cuts = []
back_m_no_cuts = []

z_g_w_cuts = [] # as a function of mass
z_g_err_w_cuts = []
z_m_w_cuts = [] # as a function of coupling constant
z_m_err_w_cuts = []

for_g_w_cuts = []
back_g_w_cuts = []
for_m_w_cuts = []
back_m_w_cuts = []


theo_for_val = []
theo_back_val = []
theo_for_err = []
theo_back_err = []

for i in range(0,n_all):
    if i < n_m and proc != "all": # constant mass, vary coupling
        # --- no cuts
        z_m_no_cuts.append([ h_probs_no_cuts[hist_names[2]][1][g + i*n_g] for g in range(0,n_g) ])
        z_m_err_no_cuts.append([ p_err_no_cuts[hist_names[2]][1][g + i*n_g] for g in range(0,n_g) ])

        for_m_no_cuts.append([ for_no_cuts[g + i*n_g + 1] for g in range(0,n_g) ])
        back_m_no_cuts.append([ back_no_cuts[g + i*n_g + 1] for g in range(0,n_g) ])

        # --- w cuts
        z_m_w_cuts.append([ h_probs_w_cuts[hist_names[2]][1][g + i*n_g] for g in range(0,n_g) ])
        z_m_err_w_cuts.append([ p_err_w_cuts[hist_names[2]][1][g + i*n_g] for g in range(0,n_g) ])

        for_m_w_cuts.append([ for_w_cuts[g + i*n_g + 1] for g in range(0,n_g) ])
        back_m_w_cuts.append([ back_w_cuts[g + i*n_g + 1] for g in range(0,n_g) ])    
        
    if i < n_g and proc != "all": # constant coupling, vary mass
        # --- no cuts
        z_g_no_cuts.append([ h_probs_no_cuts[hist_names[2]][1][i + m*n_g] for m in range(0,n_m) ])
        z_g_err_no_cuts.append([ p_err_no_cuts[hist_names[2]][1][i + m*n_g] for m in range(0,n_m) ])

        for_g_no_cuts.append([ for_no_cuts[i + m*n_g + 1] for m in range(0,n_m) ])
        back_g_no_cuts.append([ back_no_cuts[i + m*n_g + 1] for m in range(0,n_m) ])

        # --- w cuts
        z_g_w_cuts.append([ h_probs_w_cuts[hist_names[2]][1][i + m*n_g] for m in range(0,n_m) ])
        z_g_err_w_cuts.append([ p_err_w_cuts[hist_names[2]][1][i + m*n_g] for m in range(0,n_m) ])

        for_g_w_cuts.append([ for_w_cuts[i + m*n_g + 1] for m in range(0,n_m) ])
        back_g_w_cuts.append([ back_w_cuts[i + m*n_g + 1] for m in range(0,n_m) ])
        
    if (i > n_g and i > n_m) or proc == "all":
        break

# ---------- SM ----------
SM_data_probs_no_cuts = [data_no_cuts["SM"][2*i+1]/float((data_no_cuts["SM"][2*i]+data_no_cuts["SM"][2*i+1])) for i in range(0,int(len(data_no_cuts["SM"])/2))]
SM_data_probs_w_cuts = [data_w_cuts["SM"][2*i+1]/float((data_w_cuts["SM"][2*i]+data_w_cuts["SM"][2*i+1])) for i in range(0,int(len(data_w_cuts["SM"])/2))]
SM_data_err_no_cuts = [ percentage_error(data_no_cuts["SM"][2*i+1],data_no_cuts["SM"][2*i]+data_no_cuts["SM"][2*i+1]) for i in range(0,int(len(data_no_cuts["SM"])/2))]
SM_data_err_w_cuts = [ percentage_error(data_w_cuts["SM"][2*i+1],data_w_cuts["SM"][2*i]+data_w_cuts["SM"][2*i+1]) for i in range(0,int(len(data_w_cuts["SM"])/2))]

# ---------- Exclusion Limit ----------
# ---------- Printing ----------
# ---------- Plotting ----------
mark_style = ["+",".","v","^"]
cols = plt.rcParams['axes.prop_cycle'].by_key()['color']

# --- plot 1-4: SE1 and SA1 as a function of mass and coupling
if proc != "all":
    fig1, ax1 = plt.subplots(1,1)
    fig2, ax2 = plt.subplots(1,1)
    fig3, ax3 = plt.subplots(1,1)
    fig4, ax4 = plt.subplots(1,1)

    iter=0
    for i in range(0,n_all):
        if i < n_g: # loop over coupling constants
            ax1.scatter(m_values,z_g_no_cuts[i], label=str("{:.1e}".format(g_values[i-1])))
            ax1.errorbar(m_values,z_g_no_cuts[i],yerr=z_g_err_no_cuts[i],ls="")

            ax3.scatter(m_values,z_g_w_cuts[i], label=str("{:.1e}".format(g_values[i-1])))
            ax3.errorbar(m_values,z_g_w_cuts[i],yerr=z_g_err_no_cuts[i],ls="")

        if i < n_m: # loop over masses
            ax2.scatter(g_values,z_m_no_cuts[i], label=str(m_values[i-1]))
            ax2.errorbar(g_values,z_m_no_cuts[i],yerr=z_m_err_no_cuts[i],ls="")

            ax4.scatter(g_values,z_m_w_cuts[i], label=str(m_values[i-1]))
            ax4.errorbar(g_values,z_m_w_cuts[i],yerr=z_m_err_no_cuts[i],ls="")

    # ----- set title, axses, etc. 
    fig_axes = [ax1,ax2,ax3,ax4]
    figs = [fig1,fig2,fig3,fig4]

    plot_titles = ["SE1_SA1_m_no_cuts","SE1_SA1_g_no_cuts","SE1_SA1_m_w_cuts","SE1_SA1_g_w_cuts"]

    for ax in fig_axes:
        ax.set_title("S(A)=1 & S(E)=1 "+prob_text)

        if ax in [ax1,ax3]:
            x_axis = m_values
            x_len = n_m
        elif ax in [ax2,ax4]:
            x_axis = g_values
            x_len = n_g

        if ax in [ax1,ax2]:
            y_axis = SM_data_probs_no_cuts[2]
            y_err = SM_data_err_no_cuts[2]
        elif ax in [ax3,ax4]:
            y_axis = SM_data_probs_w_cuts[2]    
            y_err = SM_data_err_w_cuts[2]

        ax.plot(x_axis,[y_axis]*x_len,label="SM",color="r")
        ax.fill_between(x_axis,[y_axis-y_err]*x_len,[y_axis+y_err]*x_len,alpha=0.2,color="orange")

        # --- set x-axis
        if ax in [ax1,ax3]:
            ax.set_xlabel(m_text+" [GeV]")
                
        elif ax in [ax2,ax4]:
            ax.set_xlabel(g_text)
            ax.set_xscale("log")
            
        # --- set y-axis
        ax.set_ylabel(prob_text)
            
        # --- set legends
        if ax in [ax2,ax4]:
            figs[iter].legend(bbox_to_anchor=(0.9,1), loc="upper left", ncol=1,title=m_text+" [GeV]")
        else:
            figs[iter].legend(bbox_to_anchor=(0.9,1), loc="upper left", ncol=1,title=g_text)

        # --- set grid
        ax.grid()

        # --- save plot
        figs[iter].savefig(plot_titles[iter]+".jpg" , bbox_inches='tight', dpi=250)

        iter+=1

if bool(cos_probs_w_cuts): # if there is data_no_cuts in cos
    # --- plot 5: angular distribution data_no_cuts for mass and coupling
    # constant g
    fig5, ax5 = plt.subplots(1,1) # no cuts
    fig7, ax7 = plt.subplots(1,1) # w cuts

    # constant m
    fig6, ax6 = plt.subplots(1,1) # no cuts
    fig8, ax8 = plt.subplots(1,1) # w cuts
    
    # all electron helicities
    fig9a, ax9a = plt.subplots(1,1) # constant g
    fig9b, ax9b = plt.subplots(1,1) # constant m

    # minus electron helicities
    fig9c, ax9c = plt.subplots(1,1) # constant m

    iter = 0
    for title in titles:
        i_cut = np.where(cos_probs_no_cuts[title]==0)[0]
        x_vals= [cos_vals[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut]
        
        if title != "SM":
            a_title=np.array(title.split("_"))

            col = cols[iter]
            marker = mark_style[iter]

            g_on = 0
            m_on = 0
            if n_m > 1 and n_g > 1:
                if g_const in a_title[1] and ".0" in str(a_title[0]): # constant coupling (1*10^(-3)), vary mass
                    m_on = 1
                if m_const in a_title[0]: # constant coupling (1 GeV), vary mass
                    g_on = 1
            elif n_m > 1:
                m_on = 1
            elif n_g > 1:
                g_on = 1

            if m_on == 1:
                label_text = str(title.replace("_"+g_const,""))
            
                if proc == "ALPs" and paper_on == 1:
                    col = cols[iter+4]
                    marker = mark_style[iter]
                elif proc == "DP" and paper_on == 1:
                    col = cols[iter]
                    marker = "."
                else:
                    marker = "."
                    col = next(ax._get_lines.prop_cycler)['color']
                

                y_vals = [ cos_probs_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                y_errs = [ cos_err_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                ax5.scatter(x_vals,y_vals, label=label_text, color=col,marker=marker)
                ax5.errorbar(x_vals,y_vals,yerr=y_errs,ls="", color=col)

                y_vals = [ cos_probs_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                y_errs = [ cos_err_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                ax7.scatter(x_vals,y_vals, label=label_text, color=col,marker=marker)
                ax7.errorbar(x_vals,y_vals,yerr=y_errs,ls="", color=col)

                ax_vals = np.divide(cos_all_s[title][1],cos_all_s[title][0],where=cos_all_s[title][0]!=0.0,out=np.zeros(len(cos_all_s[title][0])))
                ax9a.scatter(cos_vals,ax_vals,label=label_text)

                iter+=1

            if g_on == 1:
                #ax6[0].scatter(cos_vals,np.round(cos_probs[title],2), label=title)
                #ax6[1].scatter(cos_vals,np.round(cos_probs_no_cuts[title],2), label=str("{:.1e}".format(float(title.replace(m_const+"_","")))))
                label_text = str("{:.1e}".format(float(title.replace(m_const+"_",""))))

                y_vals = [ cos_probs_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                y_errs = [ cos_err_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                ax6.scatter(x_vals,np.round(y_vals,2), label=label_text,marker=".")
                ax6.errorbar(x_vals,np.round(y_vals,2),yerr=y_errs,ls="")

                y_vals = [ cos_probs_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                y_errs = [ cos_err_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                ax8.scatter(x_vals,np.round(y_vals,2), label=label_text,marker=".")
                ax8.errorbar(x_vals,np.round(y_vals,2),yerr=y_errs,ls="")

                ax_vals = np.divide(cos_all_s[title][1],cos_all_s[title][0],where=cos_all_s[title][0]!=0.0,out=np.zeros(len(cos_all_s[title][0])))
                ax9b.scatter(cos_vals,ax_vals, label=label_text,marker=".")

                ax_vals = np.divide(cos_h_minus[title][1],cos_h_minus[title][0],where=cos_h_minus[title][0]!=0.0,out=np.zeros(len(cos_h_minus[title][0])))
                ax9c.scatter(cos_vals,ax_vals, label=label_text,marker=".")
  
            if proc == "all":
                if "ALP s" in title:
                    i_proc = 1
                    marker = "v"
                elif "Dark Photon" in title:
                    i_proc = 0
                    marker = "v"
                elif "ALP t/u" in title:
                    i_proc = 2
                    marker = "."

                n_list = re.findall(r'\d+',title)
                if len(n_list) == 0:
                    m_proc = 3.0
                else:
                    m_proc = float(re.findall(r'\d+',title)[0])


                y_vals = [ cos_probs_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                y_errs = [ cos_err_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
                ax5.scatter(x_vals,y_vals, label=title, color=col, marker=marker)
                ax5.errorbar(x_vals,y_vals,yerr=cos_err_no_cuts[title][0:i_cut],ls="", color=col)
                
                iter+=1

    # ----- set title, axses, etc. 
    fig_axes = [ax5,ax6,ax7,ax8,ax9a,ax9b,ax9c]
    figs = [fig5,fig6,fig7,fig8,fig9a,fig9b,fig9c]

    plot_titles = ["cos_m_no_cuts","cos_g_no_cuts","cos_m_w_cuts","cos_g_w_cuts","cos_m_all_s","cos_g_all_s","cos_h_minus"]

    iter=0
    for ax in fig_axes:   
        ax.set_ylim([0.0,1.0])
        ax.set_yticks(np.linspace(0.0,1.0,11),["0","10","20","30","40","50","60","70","80","90","100"])

        ax.grid(True)
        ax.grid(True)

        i_cut = np.where(cos_probs_no_cuts["SM"]==0)[0]
        x_vals= [cos_vals[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut]

        if ax in [ax5,ax6]:
            ax_vals = cos_probs_no_cuts["SM"]
            ax_errs = cos_err_no_cuts["SM"]

            ax.plot(cos_vals,[SM_data_probs_no_cuts[2]]*len(cos_vals),color="r",linewidth=0.8)
        elif ax in [ax7,ax8]:
            ax_vals = cos_probs_w_cuts["SM"]
            ax_errs = cos_err_w_cuts["SM"]

            ax.plot(cos_vals,[SM_data_probs_w_cuts[2]]*len(cos_vals),color="r",linewidth=0.8)
        elif ax in [ax9a,ax9b]:
            ax_list = cos_all_s
            ax_vals = np.divide(ax_list["SM"][1],ax_list["SM"][0],where=ax_list["SM"][0]!=0.0,out=np.zeros(len(ax_list["SM"][1])))
        elif ax in [ax9c]:
            ax_list = cos_h_minus
            ax_vals = np.divide(ax_list["SM"][1],ax_list["SM"][0],where=ax_list["SM"][0]!=0.0,out=np.zeros(len(ax_list["SM"][1])))

        #SM
        y_vals = [ ax_vals[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
        y_errs = [ ax_errs[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
        ax.scatter(x_vals,y_vals,label="SM",color="r",linewidth=0.8,alpha=0.8,marker="x",s=15)
        ax.errorbar(x_vals,y_vals,yerr=y_errs,ls="",color="r",linewidth=0.8,alpha=0.8)

        if proc == "all" and ax in [ax5]:
            ax.arrow(x=1,y=1,dx=0,dy=-0.3,width=.01)
            ax.annotate('Increasing \n $m_{ALP}$',xy=(0.65,0.75))

            ax.arrow(x=1,y=0,dx=0,dy=+0.3,width=.01)
            ax.annotate('Increasing \n $m_{\gamma\'}$',xy=(0.65,0.15))

        if ax in [ax5,ax7,ax9a]:
            if proc == "DP" or proc == "ALPs":
                ax.legend(loc='upper right',title=m_text+" [GeV]")
            elif proc == "ALPb" or proc == "ALPt":
                ax.legend(loc='lower left',title=m_text+" [GeV]")
            elif proc == "all":
                ax.legend(loc='lower left',title="Process")
            else:
                ax.legend(loc='upper right',bbox_to_anchor=(1.2, 1),title=m_text+" [GeV]")

        elif ax in [ax6,ax8,ax9b,ax9c]:
            if proc == "DP" or proc == "ALPs":
                ax.legend(loc='upper right',title=g_text)
            elif proc == "ALPb" or proc == "ALPt":
                ax.legend(loc='lower left',title=g_text)
            else:
                ax.legend(loc='upper right',bbox_to_anchor=(1.2, 1),title=g_text)
        else:
            ax.legend()

        ax.set_ylabel(prob_text)
        ax.set_xlabel("$cos(\\theta_{\\gamma})$")
        ax.axvline(0,color='black', linewidth=0.5, alpha=0.8) # add line at x = 0


        if proc != "all":
            if ax in [ax5,ax7]:
                #figs[iter].suptitle("$s(\\gamma) = s(e^-)$ as a function of $cos(\\theta_{\\gamma})$ and "+m_text)
                text1=str("{:.1e}".format(float(g_const))).split("e")
                text2="$"+text1[0]+"x10^"+text1[1].replace("+","").replace("0","")+"$"
                #ax.set_title(g_text+"= "+text2)
            if ax in [ax9a]:
                figs[iter].suptitle("$s(\\gamma) = +1$ as a function of $cos(\\theta_{\\gamma})$ and "+m_text)
                ax.set_title(g_text+"= "+str("{:.1e}".format(float(g_const))))
            elif ax in [ax6,ax8,ax9b,ax9c]:
                figs[iter].suptitle("$s(\\gamma) = s(e^-)$ as a function of $cos(\\theta_{\\gamma})$ and "+g_text)
                ax.set_title(m_text+"= "+m_const+" GeV")

        # --- save plot
        figs[iter].savefig(plot_titles[iter]+".svg", bbox_inches='tight')#, dpi=250)

        iter+=1

# --- plot 10: function of mass
fig10a, ax10a = plt.subplots(1,1)
fig10b, ax10b = plt.subplots(1,1)

ax10a.plot(m_values,[for_no_cuts[SM_pos][0]]*n_m,label="SM",color="r")
errMINUS=[for_no_cuts[SM_pos][0]]*n_m-for_no_cuts_errs[SM_pos]
errPLUS=[for_no_cuts[SM_pos][0]]*n_m+for_no_cuts_errs[SM_pos]
ax10a.fill_between(m_values,errMINUS,errPLUS,alpha=0.2,color="orange")

ax10b.plot(m_values,[back_no_cuts[SM_pos]]*n_m,label="SM",color="r")
errMINUS=[back_no_cuts[SM_pos][0]]*n_m-back_no_cuts_errs[SM_pos]
errPLUS=[back_no_cuts[SM_pos][0]]*n_m+back_no_cuts_errs[SM_pos]
ax10b.fill_between(m_values,errMINUS,errPLUS,alpha=0.2,color="orange")

if proc != "all":
    ax10a.plot(m_values,theo_for_val, color='orange')
    ax10b.plot(m_values,theo_back_val, color='orange')

for i in range(1,n_all):
    if i > n_g and i > n_m:
        break
    else:
        if n_m > i and i < n_g: # loop over coupling constants
            col = next(ax._get_lines.prop_cycler)['color']

            ax10a.scatter(m_values,for_g_no_cuts[i], label=str("{:.1e}".format(g_values[i])), color=col)
            ax10a.errorbar(m_values,for_g_no_cuts[i],yerr=for_no_cuts_errs[i],ls="", color=col)

            ax10b.scatter(m_values,back_g_no_cuts[i], label=str("{:.1e}".format(g_values[i])), color=col)
            ax10b.errorbar(m_values,back_g_no_cuts[i],yerr=back_no_cuts_errs[i],ls="")

if n_g > 1 and n_m > 1:
    for ax in [ax10a,ax10b]:
        # --- set x-axis
        ax.set_xlabel(m_text+" [GeV]")

        # --- set grid
        ax.grid()
        
        # --- set y-axis
        ax.set_ylabel(prob_text)

        if ax in [ax10a]:
            ax.set_title("$cos(\\theta_\\gamma) > 0$")
        elif ax in [ax10b]:
            ax.set_title("$cos(\\theta_\\gamma) < 0$")

for fig in [fig10a,fig10b]:
    fig.suptitle("S(A)=1 & S(E)=1 "+prob_text)

fig10a.legend(loc='upper left',bbox_to_anchor=(0.125, 0.8), ncol=1,title=g_text)
fig10b.legend(loc='lower left',bbox_to_anchor=(0.125, 0.1), ncol=2,title=g_text)

fig10a.savefig("for_m.jpg", bbox_inches='tight', dpi=250)
fig10b.savefig("back_m.jpg", bbox_inches='tight', dpi=250)