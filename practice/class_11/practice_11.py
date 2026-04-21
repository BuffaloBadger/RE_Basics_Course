"""Calculations for Practice Assignment 11 of REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates

# global constants available to all functions
# given
V = 500.0 # ml
tau = 0.4 # min
CA_in = 5.0E-3 # mol/ml
T_in = 60 + 273.15 # K
Cp = 1.0 # cal/ml/K
dH1 = -30.0E3 # cal/mol
k01 = 4.75E13 # /min
E1 = 25000. # cal/mol
# known
R = 1.987 # cal/mol/K
# calculated
Vdot_in = V/tau
Vdot = Vdot_in
nA_in = CA_in*Vdot_in
nZ_in = 0

# cstr model function
def cstr_model_variables(init_guess):
    # solve the cstr design equations
    soln, success, message = solve_ates(cstr_residuals, init_guess)

    # check for solver issues
    if not success:
        print('')
        print('CSTR model function issue: ' + message)
        print('')
        input('Press return to conT_inue.')

    # return the individual cstr model variables
    return soln[0], soln[1], soln[2]

# cstr residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nZ = guess[1]
    T = guess[2]

    # calculate the additional unknowns
    k1 = k01*np.exp(-E1/(R*T))
    CA = nA/Vdot
    r1 = k1*CA

    # evaluate the residuals
    epsilon_1 = nA_in - nA - V*r1
    epsilon_2 = nZ_in - nZ + V*r1
    epsilon_3 = Vdot_in*Cp*(T - T_in) + V*r1*dH1

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3])

# deliverables function
def deliverables():
    # define an initial guess for the low temperature steady state
    initial_guess = np.array([nA_in, nZ_in, T_in + 5])

    # get the cstr model variables
    nA, nZ, T_low = cstr_model_variables(initial_guess)

    # calculate the conversion
    fA_low = 100*(nA_in - nA)/nA_in

    # repeat for the high temperature steady state
    #initial_guess[2] = T_in + 100
    initial_guess[2] = T_in + 200
    nA, nZ, T_high = cstr_model_variables(initial_guess)
    fA_high = 100*(nA_in - nA)/nA_in

    # repeat for the medium temperature steady state
    initial_guess[2] = 0.5*(T_high + T_low)
    nA, nZ, T_mid = cstr_model_variables(initial_guess)
    fA_mid = 100*(nA_in - nA)/nA_in

    # tabulate, show and save the results
    results_df = pd.DataFrame({'Temperature (°C)' : np.array([T_low, T_mid, T_high]) - 273.15
                               ,'Conversion (%)' : np.array([fA_low, fA_mid, fA_high])})
    print('')
    print(results_df)
    print('')
    results_df.to_csv('practice_11_results.csv',index=False)

# execution command
if __name__ == '__main__':
    deliverables()