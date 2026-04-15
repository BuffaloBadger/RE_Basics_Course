"""Steady State Calcultions for Example 6.6.4 of REB, The Book"""

# import libraries
import numpy as np
from reb_utils import solve_ates
import pandas as pd

# given constants available to all functions 
V = 500 # cm^3
Vdot_in = 1.0 # cm^3 /s
CA_in = 0.015 # mol /cm^3
CB_in = 0.015 # mol /cm^3
T_in = 50 + 273.15 # K
Cp = 0.35*4.184 # J /g /K
rho = 0.93 # g /cm^3
dH = -20000 # J /mol
k0 = 3.24E12 # cm^3 /mol /s
E = 105000 # J /mol

# known constants available to all functions
R = 8.314 # J /mol /K

# calculated constants available to all functions
nA_in = Vdot_in*CA_in
nB_in = Vdot_in*CB_in

# cstr model function
def cstr_model_variables(init_guess):
     
	# solve the ATEs
    soln, success, message = solve_ates(cstr_residuals,init_guess)

    # check that the solution is converged
    if not(success):
        print(f"  CSTR model function error: {message}")

    # return the solution as an array
    return soln

# cstr residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nB = guess[1]
    nY = guess[2]
    nZ = guess[3]
    T = guess[4]

    # calculate the rate
    k = k0*np.exp(-E/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r = k*CA*CB

    # evaluate the residuals
    residual_1 = nA_in - nA - V*r
    residual_2 = nB_in - nB - V*r
    residual_3 = -nY + V*r
    residual_4 = -nZ + V*r
    residual_5 = -Vdot_in*rho*Cp*(T - T_in) - V*r*dH

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4
                     , residual_5])

# deliverables function
def deliverables():

    # allocate storage for the results
    ss_fA = np.zeros(3)
    ss_T = np.zeros(3)

    # set an initial guess for the low T steady state
    init_guess = np.array([nA_in, nB_in, 0.0, 0.0, 50.0 + 273.15])

    # solve the reactor design equations
    soln = cstr_model_variables(init_guess)

    # save the results
    nA = soln[0]
    ss_fA[0] = 100*(nA_in - nA)/nA_in
    ss_T[0] = soln[4] - 273.15

    # set an initial guess for the middle T steady state
    init_guess = np.array([nA_in, nB_in, 0.0, 0.0, 140.0 + 273.15])

    # solve the reactor design equations
    soln = cstr_model_variables(init_guess)

    # save the results
    nA = soln[0]
    ss_fA[1] = 100*(nA_in - nA)/nA_in
    ss_T[1] = soln[4] - 273.15

    # set an initial guess for the high T steady state
    #init_guess = np.array([nA_in, nB_in, 0.0, 0.0, 270.0 + 273.15])
    init_guess = np.array([nA_in, nB_in, 0.0, 0.0, 240.0 + 273.15])

    # solve the reactor design equations
    soln = cstr_model_variables(init_guess)

    # save the results
    nA = soln[0]
    ss_fA[2] = 100*(nA_in - nA)/nA_in
    ss_T[2] = soln[4] - 273.15

    # tabulate, show and save the results
    results_df = pd.DataFrame({'Conversion (%)':ss_fA, 'Temperature (°C)':ss_T})
    results_df = results_df.round(2)
    print(results_df)
    results_df.to_csv('example_6_6_4_steady_states.csv',index=False)

    return

# calculate the deliverables
if __name__=="__main__":
    deliverables()