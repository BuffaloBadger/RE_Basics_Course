"""Calculations for Example 6.6.3 of REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ates

# set the dpi for figures
plt.rc("savefig", dpi=300)

# constants available to all functions
# given
V = 0.5 # m^3
nA_in = 70 # mol /s
nB_in = 1500 # mol /s
Vdot_in = 40E-3 # m^3 /s
k0 = 1.2E9 # m^3 /mol /s
E = 25800*4.184 # J /mol
K0 = 4.2E-18 # m^3 /mol
dH = -22400*4.184 # J /mol
Cp_A = 412 # J /mol /K
Cp_B = 75.5 # J /mol /K
T_in = 94.7 + 273.15 # K
# known
R = 8.314 # J /mol /K

# reactor model function
def reactor_model_variables(init_guess):
	# solve the ATEs
    soln, success, message = solve_ates(residuals,init_guess)

    # check that the solution is converged
    if not(success):
        print(f"The initial temperature was NOT found: {message}")

    # return the solution
    return soln

# residuals function
def residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nB = guess[1]
    nZ = guess[2]
    T = guess[3]

    # calculate the rate
    k = k0*np.exp(-E/R/T)
    K = K0*np.exp(-dH/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    CZ = nZ/Vdot_in
    r = k*CA*CB*(1 - CZ/(K*CA*CB))

    # evaluate the residuals
    residual_1 = nA_in - nA - V*r
    residual_2 = nB_in - nB - V*r
    residual_3 = -nZ + V*r
    residual_4 = -(nA_in*Cp_A + nB_in*Cp_B)*(T-T_in) - V*r*dH

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4])

# deliverables function
def deliverables():
	# set the initial guess
    excess = 5 # K
    init_guess = np.array([nA_in, nB_in, 0.0, T_in + excess])

    # solve the reactor design equations
    solution = reactor_model_variables(init_guess)

    # calculate the conversion
    fA = 100*(nA_in - solution[0])/nA_in
    T = solution[3]

    # report the results
    print(f"Guess: {T_in + excess - 273.15} °C")
    print(f"T: {T - 273.15} °C")
    print(f"fA: {fA}%")

    return

if __name__=="__main__":
    deliverables()
