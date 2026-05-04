"""Calculations for Example 7.4.1 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
dH_1 = -101.2E3 # J/mol
k0_1 = 5.11e4 * 3600 # L/mol/h
E_1 = 74.8e3 # J/mol
T_0 = 180 + 273.15 # K
V = 1900.0 # L
CA_0 = 2.9 # mol/L
CB_0 = 3.2 # mol/L
Cp = 1.23 * 4.184 # J/g/K
rho = 1.02 * 1000.0 # g/L
t_f = 2.0 # h
# known
R = 8.314 # J/mol/K
# calculated
nA_0 = CA_0*V
nB_0 = CB_0*V

# BSTR reactor function
def BSTR_model_variables():
    # set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, nB_0, 0.0, 0.0, T_0])

    # define the stopping criterion
    f_var = 0
    f_val = t_f
     
    # solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , BSTR_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not(success):
        print('')
        print(f"BSTR function issue: {message}")
        print('')
        input('Press return to continue.')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]
    T = dep[4,:]

    # return all profiles
    return t, nA, nB, nY, nZ, T

# BSTR derivatives function
def BSTR_derivatives(ind, dep):
    # extract necessary dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    T = dep[4]

    # calculate the rate
    CA = nA/V
    CB = nB/V
    k = k0_1*np.exp(-E_1/R/T)
    r = k*CA*CB

    # evaluate the derivatives
    dnAdt = -V*r
    dnBdt = -V*r
    dnYdt = V*r
    dnZdt = V*r
    dTdt = -r*dH_1/rho/Cp

    # return the derivatives
    return [dnAdt, dnBdt, dnYdt, dnZdt, dTdt]

# deliverables function
def deliverables():
    # solve the reactor design equations
    t, nA, nB, nY, nZ, T = BSTR_model_variables()

    # calculate the other quantities of interest
    CA = nA/V
    CB = nB/V
    CY = nY/V
    CZ = nZ/V
    k = k0_1*np.exp(-E_1/R/T)
    r = k*CA*CB
    T_C = T - 273.15
    
    # display and save the graphs
    plt.figure() # concentration profiles
    plt.plot(t,CA,label='A')
    plt.plot(t,CB,label='B')
    plt.plot(t,CY, label='Y')
    plt.plot(t,CZ,linestyle=':',linewidth=4,label='Z')
    plt.xlabel("Time (h)")
    plt.ylabel("Concentration (mol L$^{-1}$)")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend()
    plt.savefig('example_7_4_1_Ci_vs_t.png')
    plt.savefig('example_7_4_1_Ci_vs_t.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex1/example_7_4_1_Ci_vs_t.png')
    plt.show(block = False)

    plt.figure() # temperature profile
    plt.plot(t,T_C)
    plt.xlabel("Time (h)")
    plt.xlim(left=0)
    plt.ylabel("Temperature (°C)")
    plt.savefig('example_7_4_1_T_vs_t.png')
    plt.savefig('example_7_4_1_T_vs_t.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex1/example_7_4_1_T_vs_t.png')
    plt.show(block = False)

    plt.figure() # rate profile
    plt.plot(t,r)
    plt.xlabel("Time (h)")
    plt.ylabel("Rate (mol L$^{-1}$ h$^{-1}$)")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.savefig('example_7_4_1_r_vs_t.png')
    plt.savefig('example_7_4_1_r_vs_t.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex1/example_7_4_1_r_vs_t.png')
    plt.show()

    return

# execution command
if __name__ == '__main__':
    deliverables()
    