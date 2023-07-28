#!/usr/bin/python
# ---------- Libaries ----------
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('tableau-colorblind10')

import sys
import os
import math
import numpy as np
from pathlib import Path
import re

from Extra.stat_extra import *
from Extra.fit_funcs import *
from Extra.DP_ALP_extra import *

# ---------- Command Line Arguments ----------
if ( len(sys.argv) < 2 ): # less than three command line arguments
    print("Not correct number of command line arguments")
    quit()

# ---------- Variables ----------
direc = sys.argv[1] # directory path
prob_text = "Percentage (%)" #"Probability"

# ---------- Mass and coupling variables
g_text = "$g_X$"
m_text = "$M_X$"

m_values = np.array([])
g_values = np.array([])

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

if (n_g > 1 ): # less than three command line arguments
    print("Too many coupling constants")
    quit()
else:
    g_const = str(g_values[0])

# mX values
m_max = max(m_values)
m_min = min(m_values)

if n_m > 1:
    m_step = m_values[1]-m_values[0] # step size
    m_plot = np.arange(m_min,m_max+1/4*m_step,step=1/4*m_step)


# ----- Sort Input parameters
m_list = np.arange(0,n_m)

if not np.all(m_values[:-1] <= m_values[1:]):
    print("m values not in order")
    m_values, m_list = (np.array(x) for x in zip(*sorted(zip(m_values,m_list), key=lambda pair:pair[0])))

# fix order of titles
titles_temp = titles.copy()
titles_temp = np.delete(titles_temp,SM_pos)

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

m_vals = []

for j in range(0,n_g*n_m):
    m_vals.append(m_values[int(j/n_g)])

#print(m_values)
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

# ----- Cos value erros ----- 
if cos_tots_w_cuts:
    for title in titles:
        cos_err_w_cuts[title] = [ percentage_error(cos_Se_Sa_w_cuts[title][c_val],cos_tots_w_cuts[title][c_val]) for c_val in range(0,len(cos_vals)) ]
        cos_err_no_cuts[title] = [ percentage_error(cos_Se_Sa_no_cuts[title][c_val],cos_tots_no_cuts[title][c_val]) for c_val in range(0,len(cos_vals)) ]
    cos_err_w_cuts["SM"] = [ percentage_error(cos_Se_Sa_w_cuts["SM"][c_val],cos_tots_w_cuts["SM"][c_val]) for c_val in range(0,len(cos_vals)) ]
    cos_err_no_cuts["SM"] = [ percentage_error(cos_Se_Sa_no_cuts["SM"][c_val],cos_tots_no_cuts["SM"][c_val]) for c_val in range(0,len(cos_vals)) ]

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

# --- plot 1: angular distribution data_no_cuts for mass and coupling
# constant g
fig1, ax1 = plt.subplots(1,1) # no cuts
fig2, ax2 = plt.subplots(1,1) # w cuts

# all electron helicities
fig3, ax3 = plt.subplots(1,1) # constant g

iter = 0
for title in titles:
    i_cut = np.where(cos_probs_no_cuts[title]==0)[0]
    x_vals= [cos_vals[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut]
    
    if title != "SM":
        a_title=np.array(title.split("_"))

        col = cols[iter]
        marker = "."
        label_text = str(title.replace("_"+g_const,""))
        
        y_vals = [ cos_probs_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
        y_errs = [ cos_err_no_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
        ax1.scatter(x_vals,y_vals, label=label_text, color=col,marker=marker)
        ax1.errorbar(x_vals,y_vals,yerr=y_errs,ls="", color=col)

        y_vals = [ cos_probs_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
        y_errs = [ cos_err_w_cuts[title][x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
        ax2.scatter(x_vals,y_vals, label=label_text, color=col,marker=marker)
        ax2.errorbar(x_vals,y_vals,yerr=y_errs,ls="", color=col)


        ax_vals = np.divide(cos_all_s[title][1],cos_all_s[title][0],where=cos_all_s[title][0]!=0.0,out=np.zeros(len(cos_all_s[title][0])))
        y_vals = [ ax_vals[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
        ax3.scatter(x_vals,y_vals,label=label_text)

        iter+=1

# ----- set title, axses, etc. 
fig_axes = [ax1,ax2,ax3]
figs = [fig1,fig2,fig3]

plot_titles = ["cos_m_no_cuts","cos_m_w_cuts","cos_m_all_s"]

iter=0
for ax in fig_axes:   
    ax.set_ylim([0.0,1.0])
    ax.set_xlim([-1.0,1.0])
    ax.set_yticks(np.linspace(0.0,1.0,11),["0","10","20","30","40","50","60","70","80","90","100"])

    ax.grid(True)
    ax.grid(True)

    i_cut = np.where(cos_probs_no_cuts["SM"]==0)[0]
    x_vals= [cos_vals[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut]

    if ax in [ax1]:
        ax_vals = cos_probs_no_cuts["SM"]
        ax_errs = cos_err_no_cuts["SM"]

        ax.plot(cos_vals,[SM_data_probs_no_cuts[2]]*len(cos_vals),color="r",linewidth=0.8)
    elif ax in [ax2]:
        ax_vals = cos_probs_w_cuts["SM"]
        ax_errs = cos_err_w_cuts["SM"]

        ax.plot(cos_vals,[SM_data_probs_w_cuts[2]]*len(cos_vals),color="r",linewidth=0.8)
    elif ax in [ax3]:
        ax_list = cos_all_s
        ax_vals = np.divide(ax_list["SM"][1],ax_list["SM"][0],where=ax_list["SM"][0]!=0.0,out=np.zeros(len(ax_list["SM"][1])))

    #SM
    y_vals = [ ax_vals[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
    y_errs = [ ax_errs[x_n] for x_n in range(0,len(cos_vals)) if x_n not in i_cut ]
    ax.scatter(x_vals,y_vals,label="SM",color="r",linewidth=0.8,alpha=0.8,marker="x",s=15)
    ax.errorbar(x_vals,y_vals,yerr=y_errs,ls="",color="r",linewidth=0.8,alpha=0.8)

    if ax in [ax1,ax2,ax3]:
        ax.legend(loc='upper right',bbox_to_anchor=(1.2, 1),title=m_text+" [GeV]")
    else:
        ax.legend()

    ax.set_ylabel(prob_text)
    ax.set_xlabel("$cos(\\theta_{\\gamma})$")
    ax.axvline(0,color='black', linewidth=0.5, alpha=0.8) # add line at x = 0


    if ax in [ax1,ax2]:
        figs[iter].suptitle("$s(\\gamma) = s(e^-)$ as a function of $cos(\\theta_{\\gamma})$ and "+m_text)
    if ax in [ax3]:
        figs[iter].suptitle("$s(\\gamma) = +1$ as a function of $cos(\\theta_{\\gamma})$ and "+m_text)
        
    text1=str("{:.1e}".format(float(g_const))).split("e")
    text2="$"+text1[0]+"x10^"+text1[1].replace("+","").replace("0","")+"$"
    ax.set_title(g_text+"= "+text2)

    # --- save plot
    figs[iter].savefig(plot_titles[iter]+".jpg", bbox_inches='tight')#, dpi=250)

    iter+=1