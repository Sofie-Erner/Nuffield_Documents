# ---------- Libaries ----------
#from matplotlib.lines import _LineStyle
import numpy as np
from numpy.ma.core import sqrt
import scipy.stats as stats
import math
import sympy
from scipy.special import binom

# ---------- Functions ----------
# for fitting results

# overall function
def find_fit(x_values,y_values,SM,cl_low,cl_up,h_min,h_max,file,degree):
    # only fit through SM if it is a part of the y_values
    zero_arr = np.where(np.array(y_values) == SM)[0]

    if len(zero_arr) != 0:
        # set weigths to force fit through SM value
        sig = [0.1]*len(x_values)
        
        # doesn't fit middle x_values well
        if zero_arr[0] == int(len(x_values)/2): # SM in the middle
            for i in range(int(len(x_values)/4)+1,int(3*len(x_values)/4)):
                if i != int(len(x_values)/2):
                    sig[i] = 10

        sig[zero_arr[0]] = 100
        
        SM_in = 1
    else:
        SM_in = 0

    if degree == 0: # default if degree not known
        if SM_in == 1:
            fit_func = np.poly1d(np.polyfit(x_values, y_values,1,w=sig),variable='h')
        elif SM_in == 0:
            fit_func = np.poly1d(np.polyfit(x_values, y_values,1),variable='g')
        SSR = np.array((fit_func(x_values) - SM)**2).sum()
        SST = ((y_values - SM)**2).sum()

        # calculate R squared
        r2 = SSR/SST
        if r2 >= 0.98: # for a good linear fit
            degree = 1
        else:
            if SM_in == 1:
                fit_func = np.poly1d(np.polyfit(x_values, y_values,2,w=sig),variable='h')
            elif SM_in == 0:
                fit_func = np.poly1d(np.polyfit(x_values, y_values,2),variable='g')

            SSR = np.array((fit_func(x_values) - SM)**2).sum()
            r2 = SSR/SST
            if r2 >= 0.98:# good quadratic fit
                degree = 2
            else: # do quartic fit
                if SM_in == 1:
                    fit_func = np.poly1d(np.polyfit(x_values, y_values,4,w=sig), variable='h')
                elif SM_in == 0:
                    fit_func = np.poly1d(np.polyfit(x_values, y_values,4), variable='g')
    elif degree in [1,2,4]:
        if SM_in == 1:
            fit_func = np.poly1d(np.polyfit(x_values, y_values,degree,w=sig),variable='h')
        elif SM_in == 0:
            fit_func = np.poly1d(np.polyfit(x_values, y_values,degree),variable='g')
    else:
        print("Unknown degree; ",degree)

    if degree == 1:
        fit_d = "Linear"
    elif degree == 2:
        fit_d = "Quadratic"
    elif degree == 4:
        fit_d = "Quartic"

    file.write(str(fit_func)+"\n")
    find_min_max(cl_low,cl_up,h_min,h_max,degree,fit_func.c,file)

    return fit_func, fit_d

def find_min_max(cl_low, cl_up, min, max, degree, array, file):
    if degree == 1: # linear
        text = "Linear"

        # Linear; y = A*h +B => h = 1/A ( y - B )
        lim_up = 1/array[0]*(cl_up - array[1]) # intersection with upper bound
        lim_low = 1/array[0]*(cl_low - array[1]) # intersection with lower bound

        if array[0] > 0: # positive gradient
            if lim_low > min: # new min
                file.write(text + ",min: " + str(lim_low)+"\n")
            if lim_up < max: # new max
                file.write(text + ",max: " + str(lim_up)+"\n")
        elif array[0] < 0: # negative gradient
            if lim_up > min: # new min
                file.write(text + ",min: " + str(lim_up)+"\n")
            if lim_low < max: # new max
                file.write(text + ",max: " + str(lim_low)+"\n")
    
    elif degree == 2: # quadratic
        text = "Quadratic"
        # Quadratic polynomial; y = A*h^2 + B*h + C => h = 1/(2A) * (-B +- sqrt(B^2 - 4*A*(C-y)))

        if cl_low == 0: # for chi squared and likelihood
            # intersections with critical chi-squared value
            like_lim_min, like_lim_max = find_quadratic(cl_up,min,max,array)

            if like_lim_max < max: # new max
                file.write(text + ",outer,max: " + str(like_lim_max)+"\n")
            if like_lim_min > min: # new min
                file.write(text + ",outer,min: " + str(like_lim_min)+"\n")
        else:
            # intersections with upper bound
            lim_up_min, lim_up_max = find_quadratic(cl_up,min,max,array)
        
            # intersections with lower bound
            lim_low_min, lim_low_max = find_quadratic(cl_low,min,max,array)

            if array[0] > 0: # upwards parabola
                # upper CL
                if lim_up_max < max: # new max
                    file.write(text + ",outer,max: " + str(lim_up_max)+"\n")
                if lim_up_min > min: # new min
                    file.write(text + ",outer,min: " + str(lim_up_min)+"\n")
            
                # lower CL
                if lim_low_max < lim_up_max: # inside global interval
                    file.write(text + ",inner,max: " + str(lim_low_max)+"\n")
                if lim_low_min > lim_up_min:
                    file.write(text + ",inner,min: " + str(lim_low_min)+"\n")

            elif array[0] < 0: # downwards parabola
                # lower CL
                if lim_low_max < max: # new max
                    file.write(text + ",outer,max: " + str(lim_low_max)+"\n")
                if lim_low_min > min: # new min
                    file.write(text + ",outer,min: " + str(lim_low_min)+"\n")
            
                # upper CL 
                if lim_up_max < lim_low_max: # inside global interval
                    file.write(text + ",inner,max: " + str(lim_up_max)+"\n")
                if lim_up_min > lim_low_min:
                    file.write(text + ",inner,min: " + str(lim_up_min)+"\n")
        
    elif degree == 4: # Quartic
        text = "Quartic"
        # Quartic polynomial; y = A*h^4 + B*h^3 + C*h^2 + D*h^1 + E
        #  1st derivative; y' = 4A*h^3 + 3B*h^2 + 2C*h + D
        # 2nd derivative; y'' = 12A*h^2 + 6B*h + 2C

        # intersections
        fit_low = array.copy()
        fit_up = array.copy()

        fit_low[4] -= cl_low 
        i_low = np.sort(np.roots(fit_low))
        fit_up[4] -= cl_up
        i_up = np.sort(np.roots(fit_up))

        # derivatives
        d1 = [4*array[0],3*array[1],2*array[2],array[3]]
        d1_roots = np.sort(np.roots(d1))
        
        n_roots = len(d1_roots)
        d2_of_roots = 12*array[0]*d1_roots**2 + 6*array[1]*d1_roots + 2*array[2]

        if array[0] > 0 and n_roots == 3:
            arr_out = i_up.copy() # array for min 1, max 2, min 3, max 4
            arr_in = i_low.copy() # array for max 1, min 2, max 3, min 4
        elif array[0] < 0 and n_roots == 3:
            arr_out = i_low.copy() # array for min 1, max 2, min 3, max 4
            arr_in = i_up.copy() # array for max 1, min 2, max 3, min 4
        else:
            print("Error in finding quartic limits: a=",a," and n_roots= ",n_roots)
            exit(1)

        for i in range(0,4):
            #root=d1_roots[i]
            #y_val = array[0]*root**4+array[1]*root**3+array[2]*root**2+array[3]*root+array[4]
            #if (cl_low > y_val or y_val > cl_up ) and y_val > 0:# and min < root < max: # if outside boundary
            if i == 0: # min 1 max 1
                quartic_print(arr_out[0],arr_in[0],min,max,1,text,array,file)
            elif i == 1: # min 2 max 2
                quartic_print(arr_in[1],arr_out[1],min,max,2,text,array,file)
            elif i == 2: # min 3 max 3
                quartic_print(arr_out[2],arr_in[2],min,max,3,text,array,file)
            elif i == 3: # min 4 max 4
                quartic_print(arr_in[3],arr_out[3],min,max,4,text,array,file)

def quartic_print(root_min,root_max,min,max,num,text,array,file):
    y_min = round(array[0]*root_min**4+array[1]*root_min**3+array[2]*root_min**2+array[3]*root_min+array[4],3)
    y_max = round(array[0]*root_max**4+array[1]*root_max**3+array[2]*root_max**2+array[3]*root_max+array[4],3)

    if min < root_min < max and not np.iscomplex(root_min) and y_min > 0: # new min
        file.write(text + ",min_"+str(num)+": " + str(root_min.real)+"\n")
        #print(text + ",min_"+str(num)+": " + str(root_min.real)+"\n")
    if min < root_max < max and not np.iscomplex(root_max)  and y_max > 0: # new max
        file.write(text + ",max_"+str(num)+": " + str(root_max.real)+"\n")
        #print(text + ",max_"+str(num)+": " + str(root_max.real)+"\n")

                
def find_quadratic(cl, min, max, array):
    lim_min = min
    lim_max = max
    # discriminant
    d = (array[1])**2 - 4*array[0]*(array[2] - cl)

    # intersections with critical value
    # d < 0 => imaginary roots
    if d > 0: # two real roots
        lim_min = 1/(2*array[0])*(-array[1] - sqrt(d))
        lim_max = 1/(2*array[0])*(-array[1] + sqrt(d))

        if lim_min > lim_max: # swap
            lim_min, lim_max = lim_max, lim_min

    elif d == 0: # one real root
        if array[0] > 0: # upwards parabola
            # everything outside confidence interval
            lim_min = 0
            lim_max = 0
        elif array[0] < 0: # downwards parabola
            # everything inside confidence interval
            lim_min = min
            lim_max = max
    
    return lim_min, lim_max
