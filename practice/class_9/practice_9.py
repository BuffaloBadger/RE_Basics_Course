"""Calculations for the Class 9 Practice Assignment in REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates

# global constants
# given
k0f = 1.93e+13 # ft3/lbmol/min
k0r = 2.86e+26 # /min
E1f = 4.64e+04 # BTU/lbmol
E1r = 8.67e+04 # BTU/lbmol
dH1 = -4.03e+04 # BTU/lbmol
nA_in = 0.15 # lbmol/min
nB_in = 3.3 # lbmol/min
Vdot_in = 1.41 # ft3/min
V = 18 # ft3
CpA = 100 # BTU/lbmol/°R
CpB = 18 # BTU/lbmol/°R
T_in = 200 + 459.67 # °R
# known
R = 1.986 # BTU /lbmol /°R

# cstr model function
def cstr_model_variables():
    # guess the cstr model variables
    guess = np.array([nA_in, nB_in, 0, T_in + 5])

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
    kf = k0f * np.exp(-E1f / R/T)
    kr = k0r * np.exp(-E1r / R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    CZ = nZ/Vdot_in
    r1 = kf*CA*CB - kr*CZ

    # evaluate the residuals
    residual_1 = nA_in - nA - V*r1
    residual_2 = nB_in - nB - V*r1
    residual_3 = -nZ + V*r1
    residual_4 = -(nA_in*CpA + nB_in*CpB)*(T-T_in) - V*r1*dH1

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4])

# deliverables function
def deliverables():
    # solve the cstr design equations
    nA, nB, nZ, T = cstr_model_variables()

    # calculate the conversion
    conversion = 100*(nA_in - nA) / nA_in

    # tabulate, show, and save the results
    data = [['T', f"{T-459.67:.0f}", '°F']
            ,['Conversion', f"{conversion:.1f}", '%']]
    results_df = pd.DataFrame(data, columns=('Item', 'Value', 'Units'))
    print("")
    print(results_df)
    results_df.to_csv("practice_9_results.csv", index=False)

# execution command
if __name__ == "__main__":
    deliverables()
