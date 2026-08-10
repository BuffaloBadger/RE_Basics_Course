"""Calculations for the Class 15 Practice Assignment from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ates
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
V = 50 # gal
T0 = 60 + 459.67 # °R
CA0 = 0.02 # lbmol /gal
Tex = 212 + 459.67 # °R
fA = 0.9
tf = 50 # min
Cp = 1 # BTU /lb
rho = 8.5 # lb /gal
dH1 = - 5150 # BTU /lbmol
k01 = 4630 # gal /lbmol /min
E1 = 6580 # BTU /lbmol
U = 37.3/60 # BTU /min /ft^2 /°R
# known
R = 1.986 # BTU /lbmol /°R
# calculated
nA0 = CA0 * V
nAf = nA0 * (1 - fA)

# global variable for the heat transfer area
g_A = float('nan')

# BSTR reactor function
def bstr_model_variables(Area):
    # make the area available to the derivatives function
    global g_A
    g_A = Area

    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nA0, 0, T0])

    # set the stopping criterion
    f_var = 0
    f_val = tf

    # solve the BSTR design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , bstr_model_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model function issue: {message}")
        print('')
        input("Press Enter to continue...")
    
    # return the bstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:]

# BSTR derivatives function
def bstr_model_derivatives(t, dep):
    # extract the dependent variables
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]

    # calculate the additional unknowns
    CA = nA / V
    k = k01 * np.exp(- E1 / (R * T))
    r = k * CA**2
    Qdot = U * g_A * (Tex - T)

    # calculate the derivatives
    dnA_dt = - r * V
    dnZ_dt = r * V
    dT_dt = (Qdot - r*V*dH1) / (rho*V*Cp)

    # return the derivatives
    return dnA_dt, dnZ_dt, dT_dt

# coupled unknown residual function
def coupled_unknown_residual(Area):
    # calculate the bstr model variables
    t, nA, nZ, T = bstr_model_variables(Area)

    # calculate the residual
    epsilon = nA[-1] - nAf

    # return the residual
    return epsilon

# deliverables function
def deliverables():
    # guess the heat transfer area
    A_guess = 1.0

    # calculate the heat transfer area
    soln, success, message = solve_ates(coupled_unknown_residual, A_guess)
    Area = soln[0]

    # check for solver issues
    if not success:
        print('')
        print(f"Coupled unknown residual function issue: {message}")
        print('')
        input("Press Enter to continue...")
    
    # solve the BSTR design equations
    t, nA, nZ, T = bstr_model_variables(Area)

    # tabulate, show, and save the results
    data = [["Area", f"{Area:.2f}", "ft^2"]
            ,["Reaction Time", f"{t[-1]:.2f}", "min"]
            ,["Final Conversion", f"{100*(nA0 - nA[-1]) / nA0:.1f}", "%"]
            ,["Final Temperature", f"{T[-1] - 459.67:.2f}", "°F"]]
    results_df = pd.DataFrame(data, columns=["Item", "Value", "Units"])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('practice_15_results.csv', index=False)

    # plot the temperature profile for discussion
    plt.figure(1)
    plt.plot(t, T - 459.67)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Temperature (°F)')
    plt.savefig('practice_15_T_vs_t.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    