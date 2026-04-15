"""Calculations for the Class 9 Learning Activity in REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates

# global constants
# given
k01 = 1.2e14 # ft3/lbmol/min
E1 = 45700 # BTU/lbmol
K01 = 6.5e-13 # ft3/lbmol
dH1 = -39700 # BTU/lbmol
Vdot_in = 0.08 # ft3/min
nA_in = 0.01 # lbmol/min
nB_in = 0.25 # lbmol/min
nZ_in = 0 # lbmol/min
V = 18 # ft3
CpA = 1000 # BTU/lbmol/degR
CpB = 180 # BTU/lbmol/degR
CpZ = 1200# BTU/lbmol/degR
T_in_values = np.array([600, 650, 700]) # degR
# known
R = 1.987 # BTU/lbmol/degR
# calculated
Vdot = Vdot_in

# allocate global storage for the current value of T_in
g_T_in = float('nan')

# cstr model function
def cstr_model_variables(T_in):
    # make T_in available to the cstr residuals function
    global g_T_in
    g_T_in = T_in

    # guess the cstr model variables
    guess = np.array([nA_in, nB_in, nZ_in, T_in + 5])

    # solve the cstr model equations
    soln, success, message = solve_ates(cstr_residuals, guess)

    # check for solver issues
    if not success:
        print(f"  CSTR model function error: {message}")
    
    # return the cstr model variables
    return soln[0], soln[1], soln[2], soln[3]

# cstr residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nB = guess[1]
    nZ = guess[2]
    T = guess[3]

    # calculate the additional unknowns
    k1 = k01 * np.exp(-E1 / (R * T))
    K1 = K01 * np.exp(-dH1 / (R * T))
    CA = nA/Vdot
    CB = nB/Vdot
    CZ = nZ/Vdot
    r1 = k1*CA*CB*(1 - CZ/(K1*CA*CB))

    # evaluate the residuals
    epsilon1 = nA_in - nA - r1*V
    epsilon2 = nB_in - nB - r1*V
    epsilon3 = nZ_in - nZ + r1*V
    epsilon4 = (nA_in*CpA + nB_in*CpB + nZ_in*CpZ)*(T - g_T_in) + r1*V*dH1

    # return the residuals
    return np.array([epsilon1, epsilon2, epsilon3, epsilon4])

# deliverables function
def deliverables():
    # allocate storage for the results
    conversion = np.ones_like(T_in_values) * float('nan')
    T_out = np.ones_like(T_in_values) * float('nan')
    # process each T_in value
    for i, T_in in enumerate(T_in_values):
        # solve the cstr design equations
        nA, nB, nZ, T_out[i] = cstr_model_variables(T_in)

        # calculate the conversion
        conversion[i] = 100*(nA_in - nA) / nA_in

    # tabulate, show, and save the results
    results_df = pd.DataFrame({'Tin (°R)': T_in_values
                , 'Tout (°R)': T_out, 'Conversion (%)': conversion})
    print("")
    print(results_df)
    results_df.to_csv("activity_9_results.csv", index=False)

# execution command
if __name__ == "__main__":
    deliverables()
