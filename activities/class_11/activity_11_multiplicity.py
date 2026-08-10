"""Calculations for the Class 11 Learning Activity for REB, The Course"""

# import libraries
import numpy as np
import matplotlib.pyplot as plt
from reb_utils import solve_ates

# set resolution for graphs
plt.rc("savefig", dpi=300)

# global constants available to all functions
# given
tau = 80. # s
P = 3. # atm
T_in = 260. + 273.15 # K
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

# allocate global storage for the current value of T
g_T = float('NaN')

# reactor function
def cstr_model_variables(init_guess):
    # solve the design equations
    soln, success, message = solve_ates(cstr_residuals, init_guess)

    # check for solver issues
    if not success:
        print('')
        print('CSTR model function issue:' + message)
        print('')
        input('Press enter to continue.')

    # return the reactor model variables
    return soln[0], soln[1], soln[2], soln[3] # nA, nB, nZ, and T

# residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nB = guess[1]
    nZ = guess[2]
    T_in = guess[3]

    # calculate the additional unknowns
    nA_in = 0.5*P*Vdot_in/Rpv/T_in
    nB_in = 0.5*P*Vdot_in/Rpv/T_in
    vDot = (nA + nB + nZ)*Rpv*g_T/P
    CA = nA/vDot
    CB = nB/vDot
    k1 = k01*np.exp(-E1/(R*g_T))
    r1 = k1*CA*CB

    # evaluate the residuals
    epsilon_1 = nA_in - nA - V*r1
    epsilon_2 = nB_in - nB - V*r1
    epsilon_3 = -nZ + V*r1
    epsilon_4 = (nA*CpA + nB*CpB + nZ*CpZ)*(g_T - T_in) + V*r1*dH1

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4])

# deliverables
def deliverables():
    # enable this function to set T
    global g_T

    # set a range of outlet temperatures
    #T_range = np.linspace(260,600, 100) + 273.15
    T_range = np.linspace(275,450, 100) + 273.15

    # set an initial guess for the cstr reactor variables
    T_in_guess = T_range[0] - 5
    nA_in_guess = 0.5*P*Vdot_in/Rpv/T_in
    nB_in_guess = 0.5*P*Vdot_in/Rpv/T_in
    init_guess = np.array([nA_in_guess, nB_in_guess, 0.0, T_in_guess])

    # allocate storage for the inlet temperatures
    T_in_range = np.ones_like(T_range) * float("NaN")

    # calculate the inlet temperature for each outlet temperature
    for i in range(0,100):
        # set the outlet temperature
        g_T = T_range[i]

        # solve the cstr design equations
        nA, nB, nZ, T_in_range[i] = cstr_model_variables(init_guess)

        # use this result as the guess for the next outlet temperature
        init_guess = np.array([nA, nB, nZ, T_in_range[i]])

    # convert to °C
    T_range = T_range - 273.15
    T_in_range = T_in_range - 273.15

    # generate a multiplicity plot
    plt.figure(1) 
    plt.plot(T_in_range, T_range,'k-')
    plt.axvline(T_in - 273.15, c='b')
    plt.xlabel("Inlet Temperature (°C)")
    plt.ylabel("Outlet Temperature (°C)")
    plt.savefig('activity_11_multiplicity_plot.png')
    plt.savefig('activity_11_multiplicity_plot.pdf')
    plt.show()
    
    return

# execution command
if __name__ == "__main__":
    deliverables()
