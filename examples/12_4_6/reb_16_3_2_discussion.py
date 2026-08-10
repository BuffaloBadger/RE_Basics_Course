"""Calculations for Example 16.3.2 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
Vpfr = 4.0 # m^3
nDotTot_0 = 1.25 # mol/s
yA_0 = 0.5
yB_0 = 0.5
T_0 = 300 # K
P = 2.5 # atm
Cp_i = 25.8 # cal/mol/K
dH = -8600 # cal/mol
k_0 = 8.12e2 # /s
E = 9500.0 # cal/mol
UA = 13.6 # cal /K /s
nDotY_0 = 0
nDotZ_0 = 0
# known
Re = 1.987 # cal/mol/K
Rw = 0.08206 * 1e-3 # m^3*atm/mol/K
# calculated
nDotA_0 = yA_0*nDotTot_0
nDotB_0 = yB_0*nDotTot_0

# make T1 available to all functions
T1 = -1.0

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nDotA = dep[0]
    nDotB = dep[1]
    nDotY = dep[2]
    nDotZ = dep[3]
    T = dep[4]

	# calculate the rate
    nDotTot = nDotA + nDotB + nDotY + nDotZ
    CA = nDotA*P/nDotTot/Rw/T
    r = k_0*np.exp(-E/Re/T)*CA

	# evaluate the derivatives
    dnDotAdV = -r
    dnDotBdV = -r
    dnDotYdV = r
    dnDotZdV = r
    dTdV = -r*dH/nDotTot/Cp_i

	# return the derivatives
    return dnDotAdV, dnDotBdV, dnDotYdV, dnDotZdV, dTdV

# pfr model
def profiles(T_1):
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nDotA_0, nDotB_0, nDotY_0, nDotZ_0, T_1])

	# define the stopping criterion
    f_var = 0
    f_val = Vpfr
     
	# solve the IVODEs
    V, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nDotA = dep[0,:]
    nDotB = dep[1,:]
    nDotY = dep[2,:]
    nDotZ = dep[3,:]
    T = dep[4,:]

    # return all profiles
    return V, nDotA, nDotB, nDotY, nDotZ, T

# perform the analysis
def perform_the_analysis():
    global T1
    # set a range of T_1 values
    T_1_range = np.linspace(T_0 + 0.05, T_0 + 100, 100)

    # allocate storage for T_2
    T_2_pfr = np.zeros(100)
    T_2_he = np.zeros(100)

    for i in range(0,100):

        # set T_1
        T1 = T_1_range[i]

        # calculate T_2 from PFR design equations
        V, nDotA, nDotB, nDotY, nDotZ, T = profiles(T1)
        T_2_pfr[i] = T[-1]

        # calculate Q for heating stream 0 
        Qdot = (nDotA_0 + nDotB_0 + nDotY_0 + nDotZ_0)*Cp_i*(T1 - T_0)

        # calculate T_3 from the heat transfer area
        T_3 = Qdot/UA/0.5 + T_0 - T_2_pfr[i] + T1

        # calculate T_2 from cooing stream 2
        T_2_he[i] = Qdot/(nDotA_0 + nDotB_0 + nDotY_0 + nDotZ_0)/Cp_i + T_3
    
    # Plot each T_2 vs. T_1
    plt.figure(1) 
    plt.plot(T_1_range - 273.15, T_2_pfr - 273.15, color = 'r', label = 'PFR')
    plt.plot(T_1_range - 273.15, T_2_he - 273.15, color = 'b', label = 'Heat Exchanger')
    plt.xlabel("$T_1$ (°C)")
    plt.ylabel("$T_2$ (°C)")
    plt.legend()
    plt.tight_layout()

    # save and show the figure
    plt.savefig('reb_16_3_2/python/multiplicity_plot.png')
    plt.show()
    # Steady states are points where curves cross

    return

if __name__=="__main__":
    perform_the_analysis()
