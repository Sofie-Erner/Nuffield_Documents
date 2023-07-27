import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors
plt.style.use('tableau-colorblind10')
import numpy as np
from pathlib import Path
import sys
import os
import math

#from scipy.special import lambertw
from scipy.stats import poisson
import scipy.special as sc

import statistics

from sympy.parsing.mathematica import parse_mathematica
from sympy.parsing.mathematica import mathematica
from sympy import var

from Extra.belle_funcs import *

# ---------- Command Line Arguments ----------
ending = ".dat"

if ( len(sys.argv) < 3 ): # less than three command line arguments
    print("Not correct number of command line arguments")
    quit()
elif ( len(sys.argv) > 3): # more than one file given
    files = sys.argv[2:] # filenames
    titles = np.array([])

    for file in files:
        indexes = [x for x, v in enumerate(file) if v == '/']

        if len(indexes) >= 1: # directory included in the file name
            title = file[(indexes[-1]+1):].replace(ending,'')
            titles = np.append(titles,title)
        else:
            titles = np.append(titles,file.replace(ending,''))

        if ending not in file:
            print("Error in file given: ", file)
            quit()
else:
    files = np.array([sys.argv[2]]) # filename
    if ending not in files[0]:
        print("Error in file given: ", files)
        quit()

direc = sys.argv[1] # directory path

# ---------- Variables ----------
SM_on = 0 # SM on or off
if len(files) > 1:
    # --- Masses & Couplings
    m_values = np.unique(np.array([ float(title.split("_")[0]) for title in titles if "SM" not in title ]))
    n_m = len(m_values)
    c_values = np.unique(np.array([ float(title.split("_")[1]) for title in titles if "SM" not in title ]))
    n_c = len(c_values)
    print(m_values)
    print(c_values)

# --- Energy & Theta Values
Emin = 1.8
Emax = 5.8
d_E = 0.1

tmin = 10
tmax = 160
d_t = 1

if len(files) == 1:
    z_fin = np.array([]) # empty array
    if "SM" in files[0]:
        SM_on = 1
else:
    z_data = {} # empty dictionary

theta_vals = np.arange(tmin,tmax+d_t,d_t)
E_vals = np.arange(Emin,Emax+d_E,d_E)
L = 20 # luminosity 20 fb^-1 50*10**3 
multi = 1 #10**4
E_ISR = 0.4

t_cols = 0
t_rows = 0
t_max = 0

# --- Photon Detection Efficiency
# approximate equations based on plot from Belle II Physics Book 
def phot_eff(m):
    if m < 6:
        return -1.451*10**(-3)*m**3 + 4.022*10**(-3)*m**2 - 2.169*10**(-3)*m + 3.166*10**(-1)
    elif m >= 6:
        return -7.383*10**(-2)*m**3 + 1.510*m**2 - 1.030*10**(+1)*m + 2.384*10**(+1)
    else:
        return 1
    

# Dark Matter Mass and CMS Energy
def E_CMS_1(x):
    E1 = 4
    E2 = 7
    s = 4*E1*E2
    return (s - x**2)/(2*np.sqrt(s))

def E_CMS_2(x):
    E1 = 4
    E2 = 7
    s = 4*E1*E2
    return (s - x)/(2*np.sqrt(s))

def M_DM_1(x):
    E1 = 4
    E2 = 7
    s = 4*E1*E2

    if type(x) == float:
        return np.sqrt(s - 2*np.sqrt(s)*x)
    else:
        mOut = np.zeros(len(x))

        for xi in range(0,len(x)):
            m2 = s - 2*np.sqrt(s)*x[xi]

            if round(m2) == 0:
                mOut[xi] = 0
            elif m2 < 0:
                mOut[xi] = -np.sqrt(-m2)
            else:
                mOut[xi] = np.sqrt(m2)

        return mOut

def M_DM_2(x):
    E1 = 4
    E2 = 7
    s = 4*E1*E2

    return (s - 2*np.sqrt(s)*x)

def M_DM_labels(y,pos):
    if y < 0:
        return ''
    elif y == 0:
        return 0
    else:
        return round(np.sqrt(y),2)


# ---------- Get Data ----------
for iter in range(0,len(files)):
    n_cols = 0
    n_rows = 0
    z_temp = np.array([])

    # -- Open File
    if len(files) == 1:
        in_file = open(os.path.join(direc, files[0]),'r')
    else:
        in_file = open(os.path.join(direc, files[iter]),'r')
    
    lines = in_file.readlines()

    for line in lines:
        line = line.replace(' \n', '')
        line = line.replace('\n', '')
        if "\t" in line:
            a_line = [float(i) for i in line.split("\t")]
        else:
            a_line = [float(i) for i in line.split(" ") if len(i) != 0]

        if len(a_line) != 0:
            z_temp = np.append(z_temp,np.array(a_line)*d_E*d_t*L)

            if n_cols == 0:
                n_cols = len(a_line)
            elif n_cols != len(a_line):
                print("Error in number of entries")
                exit
            n_rows += 1

    if len(files) == 1:
        z_fin = z_temp.reshape(n_rows,n_cols)
    else:
        if "SM" in titles[iter]:
            title = titles[iter]
            SM_on = 1
            SM_title = titles[iter]
        else:
            a_title = titles[iter].split("_") 
            title = (float(a_title[0]),float(a_title[1]))

        z_data[title] = z_temp.reshape(n_rows,n_cols)

        if np.max(z_data[title]) > t_max: # set global max for scaling
            t_max = np.max(z_data[title])
    
    if iter == 0:
        t_cols = n_cols
        t_rows = n_rows

        if t_cols != len(theta_vals) or t_rows != len(E_vals):
            print("Error: Number of rows or columns do not match given")
            quit()
    else:
        if n_cols != t_cols or n_rows != t_rows:
            print("Error: Number of rows or columns in files don't match")
            quit()

# ---------- Analysis ----------
if len(files) > 1:
    temp_titles = titles[titles != [title for title in titles if "SM" in title][0]]

    if any("SM" in title for title in titles) and "SM" not in titles[0]: # SM not first
        titles = np.append([title for title in titles if "SM" in title][0],temp_titles)

    print(titles)

# ---- Exclusion Limits
p_crit = 0.1

# --- Import fit functions
if SM_on == 1:
    """
    fitLow_file = open(r'/home/serner/Documents/Project_files/DP_ALP/MadGraph_results/SM_back/Belle/Belle_comp/LHE_results/FitLowAutoCOMB.txt','r')
    fitHigh_file = open(r'/home/serner/Documents/Project_files/DP_ALP/MadGraph_results/SM_back/Belle/Belle_comp/LHE_results/FitHighAutoCOMB.txt','r')
    in_file = [fitLow_file, fitHigh_file]
    """
    if (len(files) > 1 and "SM_all" in titles) or "SM_all" in files[0]:
        fitLow_file = open(r'/home/serner/Documents/Project_files/DP_ALP/MadGraph_results/SM_back/Belle/Belle_comp/LHE_results/FitLowAutoCOMB.txt','r')
        fitHigh_file = open(r'/home/serner/Documents/Project_files/DP_ALP/MadGraph_results/SM_back/Belle/Belle_comp/LHE_results/FitHighAutoCOMB.txt','r')
        in_file = [fitLow_file, fitHigh_file]
    elif (len(files) > 1 and "SM_RR" in titles) or "SM_RR" in files[0]:
        fit_file = open(r'/home/serner/Documents/Project_files/DP_ALP/MadGraph_results/SM_back/Belle/Belle_comp/LHE_results/FitHighAutoRR.txt','r')
        in_file = [fit_file]
        #fitLow_file = open(r'/home/serner/Documents/Project_files/DP_ALP/MadGraph_results/SM_back/Belle/Belle_comp/LHE_results/FitLowAutoRR.txt','r')
        #fitHigh_file = open(r'/home/serner/Documents/Project_files/DP_ALP/MadGraph_results/SM_back/Belle/Belle_comp/LHE_results/FitHighAutoRR.txt','r')
        #in_file = [fitLow_file, fitHigh_file]
    elif (len(files) > 1 and "SM_RL" in titles) or "SM_RL" in files[0]:
        fitLow_file = open(r'/home/serner/Documents/Project_files/DP_ALP/MadGraph_results/SM_back/Belle/Belle_comp/LHE_results/FitLowAutoCOMB.txt','r')
        fitHigh_file = open(r'/home/serner/Documents/Project_files/DP_ALP/MadGraph_results/SM_back/Belle/Belle_comp/LHE_results/FitHighAutoCOMB.txt','r')
        in_file = [fitLow_file, fitHigh_file]
    else:
        print("Error SM fit file not found")
        quit()
    #"""
    str_fit = [None]*len(in_file)
    E_low_high = [None]*len(in_file)

    for i in range(0,len(in_file)):
        lines = in_file[i].readlines()
        test1 = lines[0]
        x = var('x')
        line_fit = parse_mathematica(test1.replace("*^","*10^"))
        str_fit[i] = str(line_fit).replace("sign","np.sign").replace("Abs","np.abs").replace("sin","np.sin").replace("cos","np.cos").replace("anp.cos","np.arccos").replace("sqrt","np.sqrt")
        E_low_high[i] = R_func(theta_vals,str_fit[i])

# --- Find number of events for low and high mass cuts
if len(files) > 1 and SM_on == 1:
    # Find E_value for DM contributions
    l_temp = [ (z_data[k] != 0).sum(1) for k in z_data.keys() if "SM" not in k]
    E_DM = [ np.nonzero(l_temp[i])[0][0] for i in range(0,len(l_temp)) ]

    #print("E_DM: ",E_vals[E_DM])
    N_SM = np.zeros([2,len(E_DM)]) # number of SM events for each DM mass or coupling
    N_SM_crit = np.zeros([2,len(E_DM)])

    N_low_high = np.zeros([2,len(z_data.keys())]) # number of event for each fit
    N_crit_low_high = np.zeros(len(E_low_high)) # signal number of event for each fit

    print("---- Poisson ----")
    for n in range(0,len(E_low_high)):
        for i in range(0,len(theta_vals)):
            for j in range(0,len(E_vals)):
                if np.round(E_low_high[n][i],1) < E_vals[j]:
                    N_low_high[n] += [z_data[k][j][i] for k in z_data.keys()]
                    
                    l_temp = [ z_data[SM_title][j][i] if (E_vals[E_DM[k]]-E_ISR <= E_vals[j] and E_vals[E_DM[k]]+E_ISR >= E_vals[j] ) else 0 for k in range(0,len(E_DM)) ]
                    if len(l_temp) > 0:
                        N_SM[n] += l_temp

        # --- Use Poisson Distribution to find signal number of evetns
        N_crit_low_high[n] = sc.gammainccinv(1+int(N_low_high[n][0]),p_crit) - N_low_high[n][0]
        N_SM_crit[n] = [ sc.gammainccinv(1+int(n_SM),p_crit) - n_SM for n_SM in N_SM[n] ]

        c_crit = np.zeros([2,len(z_data.keys())])
        c_crit_indi = np.zeros([2,len(z_data.keys())])

        for n in range(1,len(z_data.keys())):
            #print(titles[n])
            m_i = float(titles[n].split("_")[0])
            c_i = float(titles[n].split("_")[1])

            if N_low_high[0][n] != 0:
                #n_SM = N_crit_low_high[0]
                #n_SM = N_SM_crit[0][n-1]
                if c_i < 1:
                    c_crit[0][n] = np.sqrt((N_crit_low_high[0]/(phot_eff(m_i)*L))/(N_low_high[0][n]/L))*c_i
                    c_crit_indi[0][n] = np.sqrt((N_SM_crit[0][n-1]/(phot_eff(m_i)*L))/(N_low_high[0][n]/L))*c_i
                else:
                    c_crit[0][n] = np.sqrt((N_crit_low_high[0]/(phot_eff(m_i)*L))/(N_low_high[0][n]/L))*(1/c_i)
                    c_crit_indi[0][n] = np.sqrt((N_SM_crit[0][n-1]/(phot_eff(m_i)*L))/(N_low_high[0][n]/L))*(1/c_i)
            
            if len(E_low_high) > 1:
                if N_low_high[1][n] != 0:
                    #n_SM = N_crit_low_high[1]
                    #n_SM = N_SM_crit[1][n-1]
                    if c_i < 1:
                        c_crit[1][n] = np.sqrt((N_crit_low_high[1]/(phot_eff(m_i)*L))/(N_low_high[1][n]/L))*c_i
                        c_crit_indi[1][n] = np.sqrt((N_SM_crit[1][n-1]/(phot_eff(m_i)*L))/(N_low_high[1][n]/L))*c_i
                    else:
                        c_crit[1][n] = np.sqrt((N_crit_low_high[1]/(phot_eff(m_i)*L))/(N_low_high[1][n]/L))*(1/c_i)
                        c_crit_indi[1][n] = np.sqrt((N_SM_crit[1][n-1]/(phot_eff(m_i)*L))/(N_low_high[1][n]/L))*(1/c_i)
   
    

    print(N_low_high)
    #print(N_crit_low_high)
    print("N_SM: ",N_SM)
    #print(N_SM_crit)
   
    print("----")
    #print([ np.sqrt((N_SM_crit[0][n]/(phot_eff(float(titles[n+1].split("_")[0]))*L))/(N_low_high[0][n+1]/L))*float(titles[n+1].split("_")[1]) for n in range(0,len(titles)-1) ])
    #print([ np.sqrt((N_SM_crit[1][n]/(phot_eff(float(titles[n+1].split("_")[0]))*L))/(N_low_high[1][n+1]/L))*float(titles[n+1].split("_")[1]) for n in range(0,len(titles)-1) ])

    #print(c_crit)

    if n_c > 1 and n_m == 1: # more than one coupling constant for one mass
        l_i = -(round(np.log10(c_crit[1]))-1)
        
        # check if each coupling gives same critical value
        if len(np.unique(np.round(c_crit[1:], l_i))) != 1:
            print("Error, did not get some critical value for mass")
            quit()

    """
    print([poisson.cdf(k=int(N_low[0]), mu=(int(N_low[0]) + N_low[n])) for n in range(0,len(N_low)) ])
    print([poisson.cdf(k=int(N_high[0]), mu=(int(N_high[0]) + N_high[n])) for n in range(0,len(N_high)) ])

    for n in range(1,len(z_data.keys())):
        print("----------",titles[n],"----------")
        print(round(np.log10(N_low[n])))
        # Low Mass
        if ( N_low[n] != 0 ):
            i = find_signal(int(N_low[0]),N_low[n],p_crit)

            print(titles[n],' LOW: final muS: ',(1+i)*N_low[n], "for multiplier: ", (1+i))

        # High Mass
        i = find_signal(int(N_high[0]),N_high[n],p_crit)
            
        print(titles[n],' HIGH: final muS: ',(1+i)*N_high[n], "for multiplier: ", (1+i))
    """

    print(c_crit)
    print(c_crit_indi)
# ---------- Plot ----------
print("---- Plot ----")
colMaps = ["Purples","Blues","Greens","Purples","Blues","Greens","Purples","Blues","Greens","Purples","Blues","Greens"]
cols = ["purple","blue","green","purple","blue","green","purple","blue","green","purple","blue","green","purple","blue","green","purple","blue","green"]


# --- Contour Map
fig1, ax1 = plt.subplots(layout='constrained')
labels = np.array([])
z_label = r'Entries/bin ($\mathcal{L}=$'+str(L)+r'$fb^{-1}$)'

if len(files) == 1:
    setMin, setMax = plot_min_max(z_fin)

    norm = colors.LogNorm(vmin=setMin, vmax=setMax)
    ima = ax1.pcolor(theta_vals,E_vals,z_fin, cmap="YlOrBr", shading='auto', norm=norm)
    cbar = fig1.colorbar(ima, label=z_label, format=ticker.FuncFormatter(fmt), ax=ax1) #[$fb^{-1}$]

    out_file = files[0].replace(ending,"")

else:
    out_file = ""

    # --- scaling
    ss = StandardScaler(with_mean=True,with_std=True) # set up scaler
    trans = MinMaxScaler(feature_range=(0,t_max))

    #z_scaled = ss.fit_transform([z_data[title].flatten('C') for title in titles])
    z_scaled = trans.fit_transform([z_data[n].flatten('C') for n in z_data.keys()])

    #z_fin = np.zeros(shape=(t_rows,t_cols))
    if n_m == 1:
        for iter in range(0,2):
            if "SM" in titles[iter]:
                title = titles[iter]
            else:
                title = (float(titles[iter].split("_")[0]),float(titles[iter].split("_")[1]))
                
            if np.max(z_data[title]) != t_max: # set global max for scaling
                z_temp = z_data[title]*(t_max/np.max(z_data[title]))
            else:
                z_temp = z_data[title]

            setMin, setMax = plot_min_max(z_temp)

            norm = colors.LogNorm(vmin=setMin, vmax=setMax)

            if "SM" in title:
                col_i = "YlOrBr"
                t_label = ""
            else:
                col_i = colMaps[iter]

                indexes = [x for x, v in enumerate(titles[iter]) if v == '_']

                if len(indexes) >= 1:
                    t_label = titles[iter][:indexes[0]]
                else:
                    t_label = titles[iter]

                labels = np.append(labels,t_label)

            ima = ax1.pcolor(theta_vals,E_vals,z_temp, cmap=col_i, shading='auto', norm=norm)

            if iter == 0: # first file
                c_label = z_label
                cbar = fig1.colorbar(ima, label=c_label, format=ticker.FuncFormatter(fmt), ax=ax1) #[$fb^{-1}$]

        out_file = "DM_" + [title for title in titles if "SM" in title][0]
        markers = [ plt.Line2D([0,0],[0,0], color=cols[iter], marker='o', linestyle='') for iter in range(0,2) if "SM" not in titles[iter] ]
        fig1.legend(markers, np.unique(labels), numpoints=1, title="Mass [GeV]", bbox_to_anchor=(1.05,1), loc='upper right')

    else:
        DM_cols = np.array([])
        for iter in range(0,len(titles)):
            #print(titles[iter])
            if "SM" in titles[iter]:
                title = titles[iter]
            else:
                title = (float(titles[iter].split("_")[0]),float(titles[iter].split("_")[1]))
                if title[0] not in [2.0,4.0,6.0,8.0]:
                    continue
                
            if np.max(z_data[title]) != t_max: # set global max for scaling
                z_temp = z_data[title]*(t_max/np.max(z_data[title]))
            else:
                z_temp = z_data[title]
                
            #z_temp = z_scaled[iter].reshape(t_rows,t_cols)
            #print([ z_temp.flatten('C')[j]/z_data[title].flatten('C')[j] for j in range(0,t_rows*t_cols) if z_data[title].flatten('C')[j] != 0 ])

            """
            if "SM" in title:
                z_fin = z_fin + z_data[title]
            else:
                z_fin = z_fin + multi*z_data[title]
            """
            setMin, setMax = plot_min_max(z_temp)

            norm = colors.LogNorm(vmin=setMin, vmax=setMax)

            if "SM" in title:
                col_i = "YlOrBr"
                t_label = "SM"
                ima = ax1.pcolor(theta_vals,E_vals,z_temp, cmap=col_i, shading='auto', norm=norm)
                cbar = fig1.colorbar(ima, label=z_label, format=ticker.FuncFormatter(fmt), ax=ax1, pad=0.15) #[$fb^{-1}$]
            """
            else:
                ax1.axhline(E_CMS_1(title[0]),color=cols[iter],ls="--",alpha=0.5)
            else:
                col_i = colMaps[iter]

                DM_cols = np.append(DM_cols,cols[iter])

                indexes = [x for x, v in enumerate(titles[iter]) if v == '_']

                if len(indexes) >= 1:
                    t_label = titles[iter][:indexes[0]]
                else:
                    t_label = titles[iter]
                labels = np.append(labels,t_label)
            """
            #ima = ax1.pcolor(theta_vals,E_vals,z_temp, cmap=col_i, shading='auto', norm=norm)


            #c_label = r'Entries/bin ($\mathcal{L}=20 \, fb^{-1}$)'
            #cbar = fig1.colorbar(ima, label=c_label, format=ticker.FuncFormatter(fmt), ax=ax1) #[$fb^{-1}$]

            """
            if iter == 0: # first file
                #out_file = out_file + title + "_"
                c_label = z_label
                cbar = fig1.colorbar(ima, label=c_label, format=ticker.FuncFormatter(fmt), ax=ax1) #[$fb^{-1}$]
            elif iter == len(titles)-1: # last file
                out_file = out_file + title
                c_label = ""
            else:
                out_file = out_file + title + "_"
                c_label = ""
            cbar = fig1.colorbar(ima, label=c_label, format=ticker.FuncFormatter(fmt), ax=ax1) #[$fb^{-1}$]
            """

        out_file = "DM_" + [title for title in titles if "SM" in title][0]
        #markers = [ plt.Line2D([0,0],[0,0], color=DM_cols[iter], marker='o', linestyle='') for iter in range(0,4) ]
        #fig1.legend(markers, np.unique(labels), numpoints=1, title="Mass [GeV]", bbox_to_anchor=(1.05,1), loc='upper right')

if SM_on == 1:
    if len(E_low_high) == 1:
        ex_cols = ['black']
        ex_labels = [""]
    else:
        ex_cols = ['r','black']
        ex_labels = ["Low-Mass","High-Mass"]

    for i in range(0,len(E_low_high)):
        ax1.plot(theta_vals, E_low_high[i],color=ex_cols[i])#,label=ex_labels[i])

ax1.set_xlabel(r"$\theta_{lab} [deg]$")
ax1.set_ylabel(r"$E_{CMS}$ [GeV]")
ax1.set_ylim(Emin,np.sqrt(4*4*7)/2)

# Exclusion Curves
# Secondary Dark Matter Mass Axis
secYax1 = ax1.secondary_yaxis('right',functions=(M_DM_1,E_CMS_1))

yticks = secYax1.yaxis.get_major_ticks()
yticks[1].set_visible(False)
"""
#secYax.set_yticks(secYax.get_yticks()[1:])
#secYax.yaxis.set_ticks(np.arange(0, M_DM_2(1.0), 10))
#secYax.set_ylim(M_DM_1(Emin),0)
#secYax.set_ylim(ax1.get_ylim())
#secYax.yaxis.set_major_formatter(ticker.FuncFormatter(M_DM_labels))
"""
secYax1.set_ylabel('$m_X$ [GeV]')

ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
secYax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax1.tick_params(which='minor')

fig1.savefig(out_file+".jpg" , bbox_inches='tight', dpi=250)

# --- Exclusion Limits
if len(files) > 1 and n_m > 1 and SM_on == 1:
    fig2, ax2 = plt.subplots(1,1)
    fig3, ax3 = plt.subplots(1,1)
    labels = np.array([])

    # ---------- Export
    out_f = open(out_file+'_exLims.txt','w')
    out_f.write(str(m_values).replace("[","").replace("]","").replace("\n",""))
    out_f.write("\n")

    out_f2 = open(out_file+'_exLims_indi.txt','w')
    out_f2.write(str(m_values).replace("[","").replace("]","").replace("\n",""))
    out_f2.write("\n")

    for n in range(0,len(E_low_high)):
        data = [(m_values[j-1], c_crit[n][j]) for j in range(1,len(c_crit[n])) if c_crit[n][j] != 0 ]
        x_vals = [ data[j][0] for j in range(0,len(data))]
        y_vals = [ data[j][1] for j in range(0,len(data))]

        ax2.plot(np.append(10**(-2),x_vals),np.append(y_vals[0],y_vals),marker="o",color=ex_cols[n],label=ex_labels[n])

        data = [(m_values[j-1], c_crit_indi[n][j]) for j in range(1,len(c_crit_indi[n])) if c_crit[n][j] != 0 ]
        x_vals = [ data[j][0] for j in range(0,len(data))]
        y_vals = [ data[j][1] for j in range(0,len(data))]
        ax3.plot(np.append(10**(-2),x_vals),np.append(y_vals[0],y_vals),marker="o",color=ex_cols[n],label=ex_labels[n])

        out_text = ex_labels[n]+": "+str(c_crit[n][1:]).replace("[","").replace("]","").replace("\n","")+"\n"
        out_f.write(out_text)

        out_text = ex_labels[n]+": "+str(c_crit_indi[n][1:]).replace("[","").replace("]","").replace("\n","")+"\n"
        out_f2.write(out_text)

    out_f.close()
    out_f2.close()

    for ax in [ax2,ax3]:
        ax.set_xscale('log')
        ax.set_xlim([10**(-2),10])

        ax.set_yscale('log')
        ax.set_ylim([10**(-4),10**(-1)])
        #ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1e'))

        ax.set_xlabel(r"$m_X$")
        ax.set_ylabel(r"$g_X$")

    if len(E_low_high) > 1:
        fig2.legend(loc='upper left', bbox_to_anchor=(0.15, 0.85))
        fig3.legend(loc='upper left', bbox_to_anchor=(0.15, 0.85))

    fig2.savefig(out_file+"_lims.jpg" , bbox_inches='tight', dpi=250)
    fig3.savefig(out_file+"_lims_indi.jpg" , bbox_inches='tight', dpi=250)
elif len(files) > 1 and SM_on == 1:
    # ---------- Export
    out_f = open(out_file+'_exLims.txt','w')
    out_f.write(str(m_values).replace("[","").replace("]","").replace("\n",""))
    out_f.write("\n")

    out_f2 = open(out_file+'_exLims_indi.txt','w')
    out_f2.write(str(m_values).replace("[","").replace("]","").replace("\n",""))
    out_f2.write("\n")

    for n in range(0,len(E_low_high)):
        out_text = ex_labels[n]+": "+str(c_crit[n][1:]).replace("[","").replace("]","").replace("\n","")+"\n"
        out_f.write(out_text)

        out_text = ex_labels[n]+": "+str(c_crit_indi[n][1:]).replace("[","").replace("]","").replace("\n","")+"\n"
        out_f2.write(out_text)

    out_f.close()
    out_f2.close()

# --- Plot w. delta E
if len(files) == 1 and SM_on == 1:
    fig4, ax4 = plt.subplots(layout='constrained')

    setMin, setMax = plot_min_max(z_fin)

    norm = colors.LogNorm(vmin=setMin, vmax=setMax)
    ima = ax4.pcolor(theta_vals,E_vals,z_fin, cmap="YlOrBr", shading='auto', norm=norm)
    cbar = fig4.colorbar(ima, label=z_label, format=ticker.FuncFormatter(fmt), ax=ax4)

    ax4.set_xlabel(r"$\theta_{lab} [deg]$")
    ax4.set_ylabel(r"$E_{CMS}$ [GeV]")
    ax4.set_ylim(Emin,np.sqrt(4*4*7)/2)

    if len(E_low_high) == 1:
        ex_cols = ['black']
        ex_labels = [""]
    else:
        ex_cols = ['r','black']
        ex_labels = ["Low-Mass","High-Mass"]

    for i in range(0,len(E_low_high)):
        ax4.plot(theta_vals, E_low_high[i],color=ex_cols[i])


    ax4.axhline(E_CMS_1(5.0),color=cols[iter],ls="--",alpha=0.5)
    ax4.fill_between(theta_vals,E_CMS_1(5.0)-E_ISR,E_CMS_1(5.0)+E_ISR,alpha=0.2,color=cols[iter])

    # Secondary Dark Matter Mass Axis
    secYax4 = ax4.secondary_yaxis('right',functions=(M_DM_1,E_CMS_1))

    yticks = secYax4.yaxis.get_major_ticks()
    yticks[1].set_visible(False)
    secYax4.set_ylabel('$m_X$ [GeV]')

    ax4.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    secYax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.tick_params(which='minor')

    fig4.savefig(out_file+"_dE.jpg" , bbox_inches='tight', dpi=250)

exit()
"""
z = np.ma.masked_where(z_data <= 0, z_data)
z_step = (np.max(z_data)-np.min(z_data))/100.
levels = np.arange(np.min(z_data), np.max(z_data)+z_step, z_step)
ima = ax1.contourf(theta_vals,E_vals,z, locator=ticker.LogLocator(base=10), cmap="YlOrBr")
cbar = fig1.colorbar(ima, format=ticker.FuncFormatter(fmt), label=r'Entries/bin [fb$^{-1}$]' )
cbar.formatter = ticker.LogFormatterExponent(base=10) # 10 is the default
cbar.update_ticks()

            if ( m_n < 6): # Low Mass
                if N_low_high[0][n] == 0:
                    print("Error in signal numbers for low mass")
                    quit()
                
                sigma_crit = (N_crit_low_high[0]*L)/phot_eff(m_n)
                    
                #print("Low: ",titles[n].split('_')[1],f" c crit: {np.sqrt(N_crit_low_high[0]/N_low_high[0][n])*c_n:.2E}")
                #print(f"{np.sqrt(N_crit_low_high[0]/(N_low_high[0][n]*phot_eff(m_n)))*c_n:.2E}")
                c_crit[n] = np.sqrt(N_crit_low_high[0]/(N_low_high[0][n]*phot_eff(m_n)))*c_n

            elif m_n >= 6: # High Mass
                if N_low_high[1][n] == 0:
                    print("Error in signal numbers for high mass")
                    quit()

                sigma_crit = (N_crit_low_high[1]*L)/phot_eff(m_n)

                #print("High: ",titles[n].split('_')[1],f" c crit: {np.sqrt(N_crit_low_high[1]/N_low_high[1][n])*c_n:.2E}")
                #print(f"{np.sqrt(N_crit_low_high[1]/(N_low_high[1][n]*phot_eff(m_n)))*c_n:.2E}")
                c_crit[n] = np.sqrt(N_crit_low_high[1]/(N_low_high[1][n]*phot_eff(m_n)))*c_n

            
            else:
                print("Error in mass: ", m_n)
                quit()
            """
