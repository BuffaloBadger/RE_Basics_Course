"""Calculations for the Class 14 Learning Activity from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
V = 1200 # L
nA_0 = 2500 # mol
nB_0 = 50000 # mol
nZ_0 = 0 # mol
T_0 = 300. # K
fA_f = 0.5
k0_1 = 6e5 # /min
E_1 = 42000 # J/mol
dH1_0 = -16500. # J/mol
CpA = 180. # J/mol/K
CpB = 70. # J/mol/K
CpZ = 225. # J/mol/K
# known
R = 8.3144 # J/mol/K
# calculated
nA_f = nA_0*(1-fA_f)

# BSTR reactor function
def bstr_model_variables():
    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nA_0, nB_0, nZ_0, T_0])

    # set the stopping criterion
    f_var = 1
    f_val = nA_f

    # solve the design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
        , bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model function issue: {message}")
        print('')
        input('Press return to continue.')

    # return the bstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:]

# BSTR derivatives function
def bstr_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nZ = dep[2]
    T = dep[3]

    # calculate the additional unknowns
    k_1 = k0_1*np.exp(-E_1/(R*T))
    CA = nA/V
    r_1 = k_1*CA
    dH1 = dH1_0 + (CpZ - CpA - CpB)*(T - 298)

    # evaluate the derivatives
    dnAdt = -V*r_1
    dnBdt = -V*r_1
    dnZdt = V*r_1
    dTdt = -V*r_1*dH1/(nA*CpA + nB*CpB + nZ*CpZ)

    # return the derivatives
    return dnAdt, dnBdt, dnZdt, dTdt

# deliverables function
def deliverables():
    # solve the design equations
    t, nA, nB, nZ, T = bstr_model_variables()

    # tabulate, show and save the deliverables
    data = [["Final time", f"{t[-1]:.1f}", "min"]
     ,["Final temperature",f"{T[-1]:.1f}", "K"]]
    results_df = pd.DataFrame(data, columns=["item", "value", "units"])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_14_results.csv',index=False)

    # calculate, plot, show, and save the instantaneous rate vs time
    r = k0_1*np.exp(-E_1/(R*T)) * nA/V
    plt.figure(1)
    plt.plot(t,r)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Instantaneous Rate (mol min$^{-1}$ L$^{-1}$)')
    plt.tight_layout()
    plt.savefig('activity_14_r_vs_t.png')
    plt.savefig('activity_14_r_vs_t.pdf')
    plt.show()

    return

# execution command
if __name__ == '__main__':
    deliverables()
    