"""Calculations for the Class 11 Learning Activity in REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates

# global constants available to all functions
# given
T_in = 260. + 273.15 # K
P = 3. # atm
tau = 80. # s
k01 = 1.26e6 #L /mol /s
E1 = 19000. # cal /mol
dH1 = -10500. # cal /mol
CpA = 42/4.184 # cal /mol /K
CpB = 122/4.184 # cal /mol /K
CpZ = 173/4.184 # cal /mol /K
# known
R = 1.987 # cal/mol
Rpv = 0.08206 # L atm /mol /K
# basis
V = 1. # L
# calculated
Vdot_in = V/tau
nA_in = 0.5*P*Vdot_in/Rpv/T_in
nB_in = 0.5*P*Vdot_in/Rpv/T_in

# cstr model function
def cstr_model_variables(guess):
    # solve the cstr design equations
    soln, success, message = solve_ates(cstr_residuals, guess)

    # check for solver issues
    if not success:
        print('')
        print('CSTR Model Function Solver Issue: ' + message)
        print('')
        input('Hit return to continue.')

    # return the cstr model variables
    return soln[0], soln[1], soln[2], soln[3] # nA, nB, nZ, T

# cstr residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nB = guess[1]
    nZ = guess[2]
    T = guess[3]

    # calculate the additional unknowns
    vDot = (nA + nB + nZ)*Rpv*T/P
    CA = nA/vDot
    CB = nB/vDot
    k1 = k01*np.exp(-E1/(R*T))
    r1 = k1*CA*CB

    # evaluate the residuals
    epsilon_1 = nA_in - nA - V*r1
    epsilon_2 = nB_in - nB - V*r1
    epsilon_3 = -nZ + V*r1
    epsilon_4 = (nA*CpA + nB*CpB + nZ*CpZ)*(T - T_in) + V*r1*dH1

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4])

# deliverables function
def deliverables():
    # solve the cstr design equations using a low temperature guess
    initial_guess = np.array([nA_in, nB_in, 0.0, 280 + 273.15])
    nA, nB, nZ, T_low = cstr_model_variables(initial_guess)
    fA_low = 100*(nA_in - nA)/nA_in

    # solve the cstr design equations using a mid temperature guess
    initial_guess[3] = 350 + 273.15
    nA, nB, nZ, T_mid = cstr_model_variables(initial_guess)
    fA_mid = 100*(nA_in - nA)/nA_in

    # solve the cstr design equations using a high temperature guess
    initial_guess[3] = 425 + 273.15
    nA, nB, nZ, T_high = cstr_model_variables(initial_guess)
    fA_high = 100*(nA_in - nA)/nA_in

    # tabulate, display, and save the results
    data = [["Low Conversion", fA_low,"%"]
            ,["Low Temperature",T_low - 273.15,"°C"]
            ,["Mid Conversion",fA_mid,"%"]
            ,["Mid Temperture",T_mid - 273.15,"°C"]
            ,["High Conversion",fA_high,"%"]
            ,["High Temperature",T_high - 273.15,"°C"]]
    results_df = pd.DataFrame(data, columns=["item","value","units"])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_11_results.csv',index=False)

# execution command
if __name__ == '__main__':
    deliverables()
