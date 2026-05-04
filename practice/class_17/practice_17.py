"""Calculations for the Class 17 Practice Assignment from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
V = 1200. # L
nA0 = 2500 # mol
nB0 = 5000 # mol
nZ0 = 0. # mol
T0 = 300. # K
fAf = 0.5
cpA = 180. # J /mol /K
cpB = 70. # J /mol /K
cpZ = 225. # J /mol /K
dH1_0 = -16500. # J /mol
k01 = 6.0e5 # /min
E1 = 42000. # J /mol
t_turn = 30 # min
# known
R = 8.3144
# calculated
nAf = fAf*nA0

# BSTR reactor function
def bstr_model_variables():
    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nA0, nB0, 0, T0])

    # set the stopping criterion
    f_var = 1
    f_val = nAf

    # solve the BSTR design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model function issue: {message}")
        print('')
        input('Press return to continue.')
    
    # return the BSTR model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:]

# BSTR derivatives function
def bstr_derivatives(t,dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nZ = dep[2]
    T = dep[3]

    # calculate the additional unknowns
    k1 = k01*np.exp(-E1/(R*T))
    CA = nA/V
    r1 = k1*CA
    dH1 = dH1_0 + (cpZ - cpA - cpB)*(T-298)

    # evaluate the derivatives
    dnAdt = -V*r1
    dnBdt = -V*r1
    dnZdt = V*r1
    dTdt = -V*r1*dH1/(nA*cpA + nB*cpB + nZ*cpZ)

    # return the derivatives
    return np.array([dnAdt, dnBdt, dnZdt, dTdt])

# deliverables function
def deliverables():
    # solve the design equations
    t, nA, nB, nZ, T = bstr_model_variables()

    # calculate the net rate
    rNet = nZ[-1]/(t[-1] + t_turn)

    # tabulate, show, and save the results
    data = [["Reaction Time", f"{t[-1]:.1f}", "min"]
            ,["Final Temperature", f"{T[-1]:.0f}", "K"]
            ,["Net Rate", f"{rNet:.1f}", "mol/min"]]
    results_df = pd.DataFrame(data,columns=("Item", "Value", "Units"))
    print('')
    print(results_df)
    print('')
    results_df.to_csv('practice_17_results.csv',index=False)

# execution command
if __name__ == '__main__':
    deliverables()
    