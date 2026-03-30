"""Calculations for Example 12.7.2 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates

# constants available to all functions
# given
yA_in = 0.1
yB_in = 0.65
yI_in = 0.25
T_in = 165 + 273.15 # K
P = 5.0 # atm
yA = 0.001
tau = 0.5 # min
k0 = 1.37E5 # m^3 /mol /min
E = 11100.0 # cal /mol
dH = -7200.0 # cal /mol
Cp_A = 7.6 # cal /mol /K
Cp_B = 8.2 # cal /mol /K
Cp_I = 4.3 # cal /mol /K
# known
Re = 1.987 # cal /mol /K
Rw = 8.206E-5 # m^3 atm /mol /K
# basis
Vdot_in = 1.0 # m^3 /min
# calculated
nA_in = yA_in*P*Vdot_in/Rw/T_in
nB_in = yB_in*P*Vdot_in/Rw/T_in
nI_in = yI_in*P*Vdot_in/Rw/T_in

# reactor model function
def reactor_model_variables(init_guess):
     
	# solve the ATEs
    soln, success, message = solve_ates(residuals,init_guess)

    # check that the solution is converged
    if not(success):
        print("")
        print(f"The solver did NOT converge: {message}")

    # return the solution
    return soln

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    V = guess[0]
    nB = guess[1]
    nI = guess[2]
    nZ = guess[3]
    T = guess[4]

    # calculate nA
    nA = yA*(nB + nI + nZ)/(1 - yA)

    # calculate the rate
    k = k0*np.exp(-E/Re/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r = k*CA*CB

    # evaluate the residuals
    residual_1 = nA_in - nA - V*r
    residual_2 = nB_in - nB - V*r
    residual_3 = nI_in - nI
    residual_4 = -nZ + V*r
    residual_5 = -(nA_in*Cp_A + nB_in*Cp_B + nI_in*Cp_I)*(T - T_in) - V*r*dH

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4, residual_5])

# perform the analysis
def deliverables():
	# set the initial guess
    init_guess = np.array([1.0, nB_in, nI_in, 0.0, T_in + 5.0])

    # solve the reactor design equations
    solution = reactor_model_variables(init_guess)

    # extract the individual results
    V = solution[0]
    T = solution[4] - 273.15

    # calculate the space time
    tau = V/Vdot_in

    # read in the results from the assignment
    results_df = pd.read_csv('example_6_6_2_results.csv')
    # drop all but the first 3 rows
    results_df = results_df.iloc[[0,1],:]

    # add these results
    results_df.loc[2] = ['tau ignoring expansion',f'{tau}','min']
    results_df.loc[3] = ['T ignoring expansion',f'{T}','°C']

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('example_6_6_2_results.csv',index=False)

    return

# calculate the deliverables
if __name__=="__main__":
    deliverables()
