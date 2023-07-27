# ---------- Libaries ----------
#from matplotlib.lines import _LineStyle
import numpy as np
from numpy.ma.core import sqrt
import scipy.stats as stats
import math
import sympy
from scipy.special import binom

# ---------- Functions ----------
# check if string has number in it
def has_num(s):
    return any(i.isdigit() for i in s)

# round up and down, and give step size
def round_y(y_min, y_max):
    n = 100
    step = 100

    # if the range is bigger
    if y_max - y_min > 10000:
        n = 1000
        step = 1000
    elif y_max - y_min > 2500:
        n = 1000
        step = 500
    
    return int(math.floor(y_min/n))*n, int(math.ceil(y_max/n))*n, step

def round_h(h_min,h_max):
    h_step = 0
    if h_max - h_min > 0.1:
        h_step = 0.05
    elif h_max - h_min > 0.05:
        h_step = 0.025
    else:
        h_step = 0.01

    return h_step

# calculate the runs test probability
def runs_prob(n1, n2, y1, y2):
    r = y1 + y2
    prob = 0

    if ( r % 2 == 0): # even number of runs
        s = r/2
        prob = 2*(binom(n1-1,s-1)*binom(n2-1,s-1))/binom(n1+n2,n1)
    elif ( r % 2 == 1): # odd number of runs
        s = (r+1)/2
        prob = (binom(n1-1,s-2)*binom(n2-1,s-1) + binom(n1-1,s-1)*binom(n2-1,s-2))/binom(n1+n2,n1)
    return prob

# creturns sum up to iter
def sum_up_to(y_values,iter):
    sum = 0
    for i in range(0,iter+1):
        sum += y_values[i]
    return sum

# swap two elements in list
def move_el(l, pos, val):
    if isinstance(l,np.ndarray):
        l=np.delete(l,np.where(l == val))
    else:
        l.remove(val)
    l = np.insert(l,pos,val)
    return l

def most_common(lst):
    return max(set(lst), key=lst.count)

# sort one variable
def sort_one_var(var_vals):
    n_var = len(var_vals)
    var_min = min(var_vals)
    var_max = max(var_vals)
    var_step = 0
    var_plot = np.zeros(n_var)

    if n_var > 1:
        var_step = var_vals[1]-var_vals[0] # step size
        var_plot = np.arange(var_min,var_max+1/4*var_step,step=1/4*var_step)

    return n_var, var_min, var_max, var_step, var_plot

# sort contour variables and the like
def sort_variables(m_values,g_values,titles):
    #n_m, m_min, m_max, m_step, m_plot = sort_one_var(m_values)
    #g_m, g_min, g_max, g_step, g_plot = sort_one_var(g_values)
    n_m = len(m_values)
    n_g = len(g_values)

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

    return n_m

# ---------- ROOT Analysis ----------
"""
# ---------- Variables 
scal = 1
scal_int = 0.005; # scaling 
hist_SM = hfile.Get("SM")
hist_SM.SetLineColor(ROOT.kBlack) # set colour

# define histograms for results 
hist_res =  ROOT.TH1F("res", "results", n_bins, x_min, x_max)
hist_cl_1 = ROOT.TH1F("cl1", "cl_1", n_bins, x_min, x_max) # lower confidence level
hist_cl_2 = ROOT.TH1F("cl2", "cl_2", n_bins, x_min, x_max) # upper confidence level
hist_cl_1.SetLineColor(ROOT.kGreen)
hist_cl_2.SetLineColor(ROOT.kRed)

# ---------- Loop over bins
j=1

while j <= n_bins:
    i = 0 # iterator
    cls = 1. # cls variable

    # histograms for TLimitDataSource
    back = ROOT.TH1F("b", "back", 1, 0, 1)
    data = ROOT.TH1F("d", "data", 1, 0, 1)
    sign = ROOT.TH1F("s", "sign", 1, 0, 1)

    back.Sumw2()
    data.Sumw2()
    sign.Sumw2()

    back.SetBinContent(1, hist_SM.GetBinContent(j)) # background and data are both SM for Monte Carlo   
    data.SetBinContent(1, hist_SM.GetBinContent(j))

    #print("Bin content: ", hist_SM.GetBinContent(j)) 

    # ----- Find CLs ------
    while cls > 0.05:
        i+=1
		# ----- Scaling
        if ( i == 1):
            scal = scal_int # set first scaling
        else:
            scal = ( i*scal_int)/ ((i-1)*scal_int) # further scaling is next divided by current		

        sign.SetBinContent(1, i*scal_int*hist_SM.GetBinContent(j)) # set signal to SM scaled

        mydatasource = ROOT.TLimitDataSource()
        mydatasource.AddChannel(sign, back, data) # add signal, background, and data
        myconfidence = ROOT.TLimit.ComputeLimit(mydatasource, 50000, True) # find confidence level

        cls = myconfidence.CLs()
        #print("For scaling: ", i*scal_int,": CLs= ",cls,", signal= ",sign.GetBinContent(1),", SM= ",hist_SM.GetBinContent(1))

    #print("For Bin: ",j," Scaling= ",i*scal_int,", CLs= ",cls,", SM: ",hist_SM.GetBinContent(j))
    hist_res.SetBinContent(j, sign.GetBinContent(1)) # results = signal
    hist_cl_1.SetBinContent(j, hist_SM.GetBinContent(j) - hist_res.GetBinContent(j)) # SM - confidence level = lower bound
    hist_cl_2.SetBinContent(j, hist_SM.GetBinContent(j) + hist_res.GetBinContent(j)) # SM + confidence level = upper bound

    if hist_cl_2.GetBinContent(j) > y_max:
        y_max = hist_cl_2.GetBinContent(j) # update max y value
    elif hist_cl_1.GetBinContent(j) < y_min:
        y_min = hist_cl_1.GetBinContent(j)

    back.Delete()
    data.Delete()
    sign.Delete()

    j+=1

# ---------- BSM ----------
hist_BSM_1 = hfile.Get(titles[1])
hist_BSM_1.SetLineColor(ROOT.kBlue)

hist_BSM_2 = hfile.Get(titles[2])
hist_BSM_2.SetLineColor(ROOT.kCyan)
"""

# ----- ROOT Results
"""
c1 = ROOT.TCanvas("c1", "c1", 1200, 800)
c1.DrawFrame(x_min,y_min,x_max,y_max) # draw frame so that all histograms fit

hs = ROOT.THStack("hs", "SM w CLs")
hs.Add(hist_SM, "nostack")
hs.Add(hist_cl_1, "nostack")
hs.Add(hist_cl_2, "nostack")
hs.Add(hist_BSM_1, "nostack")
hs.Add(hist_BSM_2, "nostack")

hfile.WriteObject(hs, "SM w CLs")

hist_cl_2.GetXaxis().SetTitle(spec_x)
hist_cl_2.GetYaxis().SetTitle(spec_y)

hist_cl_2.Draw("C SAME")
hist_cl_1.Draw("C SAME")
hist_SM.Draw("C SAME")
hist_BSM_1.Draw("C SAME")
hist_BSM_2.Draw("C SAME")

legend = ROOT.TLegend(0.4,0.2,0.6,0.4)
legend.AddEntry(hist_cl_1  ,"SM - signal")
legend.AddEntry(hist_cl_2,"SM + signal")
legend.AddEntry(hist_SM,"SM")
legend.AddEntry(hist_BSM_1,titles[1])
legend.AddEntry(hist_BSM_2,titles[2])
legend.Draw("SAME")

c1.Update()
c1.cd()
c1.Modified()

root_plot_title="Confidence Levels;"+spec_x+";"+spec_y
c1.SetTitle(root_plot_title)

plot_file_root = plot_file
plot_file_root = plot_file_root.replace("split_hist", "ROOT") # remove start
c1.SaveAs(plot_file_root)
plot_file_root = plot_file_root.replace("jpg", "root") # remove start
c1.SaveAs(plot_file_root)
c1.Close()
"""

# plot 6 ROOT version
"""
root_plot_title="Fit for Exclusion Limits;"+spec_x+";"+spec_y

c2 = ROOT.TCanvas("c2", "c2", 1200, 800)
c2.DrawFrame(min(h_values),y_min,max(h_values),y_max) # draw frame so that all histograms fit
c2.SetTitle(root_plot_title)

hist_SM1 = ROOT.TH1F("SM1", "SM1", n_h, min(h_values), max(h_values)) # lower confidence level
hist_SM2 = ROOT.TH1F("SM2", "SM2", n_h, min(h_values), max(h_values)) # upper confidence level

for j in range(0,n_h):
    hist_SM1.Fill(h_values[j],y_values["SM"][int(min(cross_bins))] - sqrt(chi_crit_1*y_values["SM"][int(min(cross_bins))]))
    hist_SM2.Fill(h_values[j],y_values["SM"][int(min(cross_bins))] + sqrt(chi_crit_1*y_values["SM"][int(min(cross_bins))]))

#SM_val = [y_values["SM"][int(min(cross_bins))]]*n_h
#gr1 = ROOT.TGraph(n_h,h_values,y_values["SM"][int(min(cross_bins))])
#gr1.SetTitle("Fit for Exclusion Limits")
#val = int(100*sqrt(chi_crit_1*y_values["SM"][int(min(cross_bins))]) + 5)
#gr1.SetLineWidth(val)
#gr1.SetFillStyle(3005)
#gr1.SetLineColor(2)
#gr1.Draw("*")
hist_SM1.Draw("C")
hist_SM2.Draw("C SAME")

legend = ROOT.TLegend(0.4,0.2,0.6,0.4)
legend.AddEntry(hist_SM1,"95% Confidence Interval")
legend.Draw("SAME")

c2.Update()
c2.cd()
c2.Modified()

plot_file_root = plot_file
plot_file_root = plot_file_root.replace("split_hist", "ROOT_fit") # remove start
c2.SaveAs(plot_file_root)
plot_file_root = plot_file_root.replace("jpg", "root") # remove start
c2.SaveAs(plot_file_root)
c2.Close()
"""

"""
sfig61.scatter(h_values,i_y_values[x_values[int(min(cross_bins))]], label="BSM", color='r', linewidth=0.5)
sfig61.plot(h_values, h_fit_11(h_values), '--', color='b', label="Linear fit")
sfig61.plot(h_values, h_fit_12(h_values), '--', color='g', label="Quadratic fit")
sfig61.fill_between(h_values, [chi_cl[0]]*n_h, [chi_cl[1]]*n_h, label="95% Confidence Interval", color='orange', alpha=0.2)
#sfig61.plot(h_values, [chi_cl[0]]*n_h, label="Low CL", color='orange')
#sfig61.plot(h_values, [chi_cl[1]]*n_h, label="High CL", color='orange')
sfig61.plot(h_values,[y_values["SM"][int(min(cross_bins))]]*n_h, label="SM value", color='k')

sfig62.scatter(h_values,i_y_values[x_values[int(max(cross_bins))]], label="BSM", color='r', linewidth=0.5)
sfig62.plot(h_values, h_fit_21(h_values), '--', color='b', label="Linear fit")
sfig62.plot(h_values, h_fit_22(h_values), '--', color='g', label="Quadratic fit")
sfig62.fill_between(h_values, [chi_cl[2]]*n_h, [chi_cl[3]]*n_h, label="95% Confidence Interval", color='orange', alpha=0.2)
#sfig62.plot(h_values, [chi_cl[2]]*n_h, label="Low CL", color='orange')
#sfig62.plot(h_values, [chi_cl[3]]*n_h, label="High CL", color='orange')
sfig62.plot(h_values,[y_values["SM"][int(max(cross_bins))]]*n_h, label="SM value", color='k')
"""

# for second bin
"""
fit_file.write("x value: " + str(x_values[int(max(cross_bins))]))
fit_file.write(str(h_fit_21)+"\n")
fit_file.write(str(h_fit_22)+"\n")
if r2_21 >= 0.95:
    lim_lim = 1/h_fit_21.c[0]*(chi_cl[2] - h_fit_21.c[1])
    fit_file.write("Linear, lower: "+str(lim_lim)+"\n")

    lim_lim = 1/h_fit_21.c[0]*(chi_cl[3] - h_fit_21.c[1])
    fit_file.write("Linear, upper: "+str(lim_lim)+"\n")
else:
    # Quadratic polynomial; y = A*h^2 + B*h + C => h = 1/(2A) * (-B +- sqrt(B^2 - 4*A*(C-y)))
    lim_lim = 1/(2*h_fit_22.c[0])*(-h_fit_22.c[1] - sqrt((h_fit_22.c[1])**2 - 4*h_fit_22.c[0]*(h_fit_22.c[2] - chi_cl[2])))
    fit_file.write("Quadratic, lower, 1st root: "+str(lim_lim)+"\n")
    lim_lim = 1/(2*h_fit_22.c[0])*(-h_fit_22.c[1] + sqrt((h_fit_22.c[1])**2 - 4*h_fit_22.c[0]*(h_fit_22.c[2] - chi_cl[2])))
    fit_file.write("Quadratic, lower, 2nd root: "+str(lim_lim)+"\n")
    
    lim_lim = 1/(2*h_fit_22.c[0])*(-h_fit_22.c[1] - sqrt((h_fit_22.c[1])**2 - 4*h_fit_22.c[0]*(h_fit_22.c[2] - chi_cl[3])))
    fit_file.write("Quadratic, upper, 1st root: "+str(lim_lim)+"\n")
    lim_lim = 1/(2*h_fit_22.c[0])*(-h_fit_22.c[1] + sqrt((h_fit_22.c[1])**2 - 4*h_fit_22.c[0]*(h_fit_22.c[2] - chi_cl[3])))
    fit_file.write("Quadratic, upper, 2nd root: "+str(lim_lim)+"\n")
"""


# log likelihood functions
"""
del_log_like_1a = [0.0]*(n_all) # difference to SM Extended Poisson
del_log_like_1b = [0.0]*(n_all) # difference to SM Extended Poisson, w. bin probability
del_log_like_2a = [0.0]*(n_all) # difference to SM Extended Multinomial
del_log_like_2b = [0.0]*(n_all) # difference to SM Extended Multinomial, w. bin probability
del_log_like_3a = [0.0]*(n_all) # difference to SM Multinomial
del_log_like_3b = [0.0]*(n_all) # difference to SM Multinomial, w. bin probability
"""
"""
    del_log_like_1a[i] = y_sum[0]*math.log(y_sum[0]/y_sum[i]) + 2*(y_sum[i] - y_sum[0]) + (y_values[titles[0]]*np.log(y_values[titles[0]]/y_values[titles[i]])).sum()
    del_log_like_1b[i] = y_sum[0]*math.log(y_sum[0]/y_sum[i]) + 2*(y_sum[i] - y_sum[0]) + (y_values[titles[0]]*np.log(y_values[titles[0]]/(y_sum[0]*bin_prob[titles[i]]))).sum()
    del_log_like_2a[i] =  (y_sum[i] - y_sum[0]) + (y_values[titles[0]]*np.log(y_values[titles[0]]/y_values[titles[i]])).sum()
    del_log_like_2a[i] =  (y_sum[i] - y_sum[0]) + (y_values[titles[0]]*np.log(y_values[titles[0]]/(y_sum[0]*bin_prob[titles[i]]))).sum()
    del_log_like_3a[i] = (y_values[titles[0]]*np.log(y_values[titles[0]]/y_values[titles[i]])).sum()
    del_log_like_3b[i] = (y_values[titles[0]]*np.log(y_values[titles[0]]/(bin_prob[titles[i]]*y_sum[0]))).sum()
    """
"""
del_log_like_1a = move_el(del_log_like_1a,zero_pos,del_log_like_1a[0])
del_log_like_1b = move_el(del_log_like_1b,zero_pos,del_log_like_1b[0])
del_log_like_2a = move_el(del_log_like_2a,zero_pos,del_log_like_2a[0])
del_log_like_2b = move_el(del_log_like_2b,zero_pos,del_log_like_2b[0])
del_log_like_3a = move_el(del_log_like_3a,zero_pos,del_log_like_3a[0])
del_log_like_3b = move_el(del_log_like_3b,zero_pos,del_log_like_3b[0])
"""

# ---------- Pearson's R & R^2 
"""
sums = [0.0]*5 # s_h, s_h2, s_y, s_y2, s_hy
sums[0] = h_values.sum() #s_hy = sum over h values
sums[1] = np.array(h_values**2).sum() #s_h2 = sum over h**2 values
sums[2] = np.array(i_y_values[max_diff_x_val[0]]).sum() #s_y = sum over y values
sums[3] = (np.array(i_y_values[max_diff_x_val[0]])**2).sum() #s_y2 = sum over y**2 values
sums[4] = np.array(h_values_wz*i_y_values[max_diff_x_val[0]]).sum() #s_hy = sum over h*y values
pears_r = (n_h*sums[4] -  sums[0]*sums[2])/((n_h*sums[1] - ( sums[0])**2)*(n_h*sums[3] - (sums[2])**2))
"""

# Runs
"""
M = [0.0]*n_h # number of positive differences
N = [0.0]*n_h # number of negative differences
R = [0.0]*n_h # number of runs
no_diff = [0.0]*n_h # number of no difference
runs_M = {} # positive runs
runs_N = {} # negative runs
runs_tot = {} # runs
prob_runs = [0.0]*n_h # probability for runs
#runs_mean = [0.0]*n_h # mean for runs
#runs_var = [0.0]*n_h # standard deviation for runs
runs_Z = [0.0]*n_h # Z values for runs
"""
# Runs
"""
        r_count = 1 # placeholder for number of runs
        

        for d in range(0,len(diff_values[title])+1): # loop over differences
            on_off = 0 # for adding to lists

            if d != len(diff_values[title]):
                diff_sum[iter] += diff_values[title][d] # add to sum of differences

            if d == len(diff_values[title]): # for the last element
                on_off = 1 # add to run lists
            elif d != 0: # not for the first element
                if np.sign(diff_values[title][d]) == np.sign(diff_values[title][d-1]): # if this and the previous have the same sign
                    r_count += 1 # increase R
                else:
                    on_off = 1     
            
            if (on_off == 1):
                runs_tot[title] = np.append(runs_tot[title],r_count) # for total list

                if np.sign(diff_values[title][d-1]) == 1.0: # for positive sign
                    runs_M[title] = np.append(runs_M[title],r_count)
                elif np.sign(diff_values[title][d-1]) == -1.0: # for negative sign
                    runs_N[title] = np.append(runs_N[title],r_count)

                r_count = 1 # reset R

        prob_runs[iter] = runs_prob(M[iter],N[iter],len(runs_M[title]),len(runs_N[title])) # find probability of number of runs

        #runs_mean[iter] = (2*M[iter]*N[iter])/(M[iter]+N[iter]+1)
        R_mean = (2*M[iter]*N[iter])/(M[iter]+N[iter]+1)
        #runs_var[iter] = (2*M[iter]*N[iter]*(2*M[iter]*N[iter]-M[iter]-N[iter]))/(((M[iter]+N[iter])**2)*(M[iter]+N[iter]-1))
        R_var = 2*M[iter]*N[iter]*(2*M[iter]*N[iter]-M[iter]-N[iter])
        R_var = R_var/(((M[iter]+N[iter])**2)*(M[iter]+N[iter]-1))
        #runs_Z[iter] = (len(runs_tot[title]) - R_mean)/R_var
    
        R[iter] = len(runs_tot[title]) # add number of runs to list
        """


# writing exclusion limits to file
"""
if array[0] > 0 and n_roots == 3:
            # outer limits
            if min < i_up[3] < max and not np.iscomplex(i_up[3]): # new max
                file.write(text + ",max_4: " + str(i_up[3].real)+"\n")
            if min < i_up[0] < max and not np.iscomplex(i_up[0]): # new min
                file.write(text + ",min_1: " + str(i_up[0].real)+"\n")

            for i in range(0,len(d1_roots)): # loop over extrema
                root=d1_roots[i]
                y_val = array[0]*root**4+array[1]*root**3+array[2]*root**2+array[3]*root+array[4]

                if (cl_low > y_val or y_val > cl_up ) and y_val > 0:# and min < root < max: # if outside boundary
                    if i == 0 and not np.iscomplex(i_low[2]) and not np.iscomplex(i_low[1]): # first root
                        # has max 1 and min 2
                        if  min < i_low[0] < max and not np.iscomplex(i_low[0]): # new max
                            file.write(text + ",max_1: " + str(i_low[0].real)+"\n")
                        if  min < i_low[1] < max and not np.iscomplex(i_low[1]): # new min
                            file.write(text + ",min_2: " + str(i_low[1].real)+"\n")
                    elif i == 1 and not np.iscomplex(i_up[1]) and not np.iscomplex(i_up[2]): # middle root
                        # has max 2 and min 3
                        if  min < i_up[1] < max and not np.iscomplex(i_up[1]): # new max
                            file.write(text + ",max_2: " + str(i_up[1].real)+"\n")
                        if min < i_up[2] < max and not np.iscomplex(i_up[2]): # new min
                            file.write(text + ",min_3: " + str(i_up[2].real)+"\n")
                    elif i == 2 and not np.iscomplex(i_low[2]) and not np.iscomplex(i_low[3]): # last root
                        # has max 3 and min 4
                        if  min < i_low[2] < max and not np.iscomplex(i_low[2]): # new max
                            file.write(text + ",max_3: " + str(i_low[2].real)+"\n")
                        if min < i_low[3] < max and not np.iscomplex(i_low[3]) and i_up[3] < max: # new min
                            file.write(text + ",min_4: " + str(i_low[3].real)+"\n")

        elif array[0] < 0 and n_roots == 3:
            # outer limits
            if min < i_low[3] and i_low[3] < max and not np.iscomplex(i_low[3]): # new max
                file.write(text + ",max_4: " + str(i_low[3].real)+"\n")
            if min < i_low[0] and i_low[0] < max and not np.iscomplex(i_low[0]): # new min
                file.write(text + ",min_1: " + str(i_low[0].real)+"\n")

            for i in range(0,len(d1_roots)): # loop over extrema
                root=d1_roots[i]
                y_val = array[0]*root**4+array[1]*root**3+array[2]*root**2+array[3]*root+array[4]

                if (cl_low > y_val or y_val > cl_up ) and y_val > 0:# and min < root < max: # if outside boundary
                    if i == 0 and not np.iscomplex(i_up[2]) and not np.iscomplex(i_up[1]): # first root
                        # has max 1 and min 2
                        if  min < i_up[0] < max and not np.iscomplex(i_up[0]): # new max
                            file.write(text + ",max_1: " + str(i_up[0].real)+"\n")
                        if  min < i_up[1] < max and not np.iscomplex(i_up[1]): # new min
                            file.write(text + ",min_2: " + str(i_up[1].real)+"\n")
                    elif i == 1 and not np.iscomplex(i_low[1]) and not np.iscomplex(i_low[2]): # middle root
                        # has max 2 and min 3
                        if  min < i_low[1] < max and not np.iscomplex(i_low[1]): # new max
                            file.write(text + ",max_2: " + str(i_low[1].real)+"\n")
                        if min < i_low[2] < max and not np.iscomplex(i_low[2]): # new min
                            file.write(text + ",min_3: " + str(i_low[2].real)+"\n")
                    elif i == 2 and not np.iscomplex(i_up[2]) and not np.iscomplex(i_up[3]): # last root
                        # has max 3 and min 4
                        if  min < i_up[2] < max and not np.iscomplex(i_up[2]): # new max
                            file.write(text + ",max_3: " + str(i_up[2].real)+"\n")
                        if min < i_up[3] < max and not np.iscomplex(i_up[3]) and i_up[3] < max: # new min
                            file.write(text + ",min_4: " + str(i_up[3].real)+"\n")
"""