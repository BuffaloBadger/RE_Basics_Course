"""Calculations for Practice Assignment 11 of REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ates

# set resolution for the graph
plt.rc('savefig', dpi=300)

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

# define a global variable to make T available to all functions
g_T = float('NaN')

# cstr model function
def cstr_model_variables(init_guess):
    # solve the cstr design equations
    soln, success, message = solve_ates(cstr_residuals, init_guess)

    # check for solver issues
    if not success:
        print('')
        print('CSTR model function issue: ' + message)
        print('')
        input('Press return to continue.')

    # return the individual cstr model variables
    return soln[0], soln[1], soln[2] # nA, nZ, T_in

# cstr residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nZ = guess[1]
    T_in = guess[2]

    # calculate the additional unknowns
    k1 = k01*np.exp(-E1/(R*g_T))
    CA = nA/Vdot
    r1 = k1*CA

    # evaluate the residuals
    epsilon_1 = nA_in - nA - V*r1
    epsilon_2 = nZ_in - nZ + V*r1
    epsilon_3 = Vdot_in*Cp*(g_T - T_in) + V*r1*dH1

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3])

# deliverables function
def deliverables():
    # define a range of outlet temperatures
    #T_range = np.linspace(25, 400, 100) + 273.15
    T_range = np.linspace(25, 250, 100) + 273.15

    # allocate storage for the corresponding inlet temperatures
    T_in_range = np.ones_like(T_range) * float('NaN')

    # define an initial guess for the first outlet temperature in the range
    init_guess = np.array([nA_in, nZ_in, T_range[0] - 5])

    # loop through all of the outlet temperatures
    for i in range(100):
        # make T available to the residuals function
        global g_T
        g_T = T_range[i]

        # solve the design equations for the first outlet temperature
        nA, nZ, T_in_range[i] = cstr_model_variables(init_guess)

        # save the guess for the next outlet temperature
        init_guess = np.array([nA, nZ, T_in_range[0]])

    # generate, show, and save the multiplicity plot
    plt.figure(1)
    plt.plot(T_in_range - 273.15, T_range - 273.15, color = 'k')
    plt.axvline(T_in - 273.15, color = 'b')
    plt.xlabel('Inlet Temperature (°C)')
    plt.ylabel('Outlet Temperature (°C)')
    plt.savefig('practice_11_multiplicity_plot.pdf')
    plt.savefig('practice_11_multiplicity_plot.png')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()