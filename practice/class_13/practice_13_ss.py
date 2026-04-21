"""Calculations for Practice Assignment 13 of REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates

# global constants available to all functions
# given
V = 17.7 # ft^3
CA_in = 0.125 # lbmol/ft^3
Vdot_in = 0.15 # ft^3/min
T_in = 300 + 458.67 # °R
Cp = 73.1 # BTU/ft^3/°F
k01 = 3.04e8 # ft^3/lbmol/min
E1 = 29100 # BTU/lbmol
dH1 = 32800 # BTU/lbmol
rho_ex = 81.2 # lb/ft^3
Cp_ex = 1.06 # BTU/lb/°F
mDotex = 15.0 # lb/min
Tex_in = 410 + 458.67 # °R
Vex = 7.0 # ft^3
A = 21.5 # ft^2
U = 1.71 # BTU/ft^2/min/°F
# known
R = 1.986 # BTU/lbmol/°R
# calculated constants
Vdot = Vdot_in
nA_in = CA_in*Vdot_in

# cstr reactor model
def cstr_model_variables():
    # guess the model variables
    init_guess = np.array([nA_in, 0.0, T_in + 5, Tex_in - 5])

    # solve the cstr design equations
    soln, success, message = solve_ates(cstr_residuals, init_guess)

    # check for solver issues
    if not success:
        print('')
        print(f"CSTR model function issue: {message}")
        print('')
        input('Press return to continue.')
    
    # return the cstr model variables
    return soln[0], soln[1], soln[2], soln[3] # nA, nZ, T, Tex

# cstr residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nZ = guess[1]
    T = guess[2]
    Tex = guess[3]

    # calculate the additional unknowns
    k1 = k01*np.exp(-E1/(R*T))
    CA = nA/Vdot
    r1 = k1*CA**2
    Q = U*A*(Tex - T)

    # evaluate the residuals
    epsilon_1 = nA_in - nA - r1*V
    epsilon_2 = -nZ + r1*V
    epsilon_3 = Vdot_in*Cp*(T - T_in) + r1*V*dH1 - Q
    epsilon_4 = -Q - mDotex*Cp_ex*(Tex - Tex_in)

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4])

# deliverables function
def deliverables():
    # get the cstr model variables
    nA, nZ, T, Tex = cstr_model_variables()

    # calculate the conversion
    fA = 100*(nA_in - nA)/nA_in

    # tabulate, show, and save the results
    data =[["conversion", fA, "%"]
           ,["T", T - 458.67, "°F"]
           ,["Tex", Tex - 458.67, "°F"]]
    results_df = pd.DataFrame(data, columns=["item", "value", "units"])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('practice_13_ss_results.csv',index=False)

#execution command
if __name__=='__main__':
    deliverables()