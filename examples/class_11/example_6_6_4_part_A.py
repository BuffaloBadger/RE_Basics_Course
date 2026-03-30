"""Multiplicity plot for Example 6.6.4 of REB, The Book"""

# import libraries
import numpy as np
from reb_utils import solve_ates
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi=300)

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

# define a global variable to make T available to all functions
T = float('NaN')

# reactor model function
def reactor_model_variables(init_guess):
     
	# solve the ATEs
    soln, success, message = solve_ates(residuals,init_guess)

    # check that the solution is converged
    if not(success):
        print(f"solve_ates was unsuccessful: {message}")

    # return the solution as an array
    return soln

# residuals function
def residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nB = guess[1]
    nY = guess[2]
    nZ = guess[3]
    T_in = guess[4]

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

    # enable this function to set T
    global T
    
    # set a range of outlet temperatures
    T_range = np.linspace(-25.0, 325.0, 350) + 273.15

	# set an initial guess for the outlet molar flows and the inlet temperature
    init_guess = np.array([nA_in, nB_in, 0.0, 0.0 ,T_range[0] - 1.0])

    # allocate storage for the inlet temperatures
    Tin_range = np.zeros(350)

    # calculate the inlet temperature for each outlet temperature
    for i in range(0,350):

        # set the outlet temperature
        T = T_range[i]

        # solve the design equations
        solution = reactor_model_variables(init_guess)

        # save the inlet temperature
        Tin_range[i] = solution[4]

        # use the result as the initial guess for the next outlet T
        init_guess = solution

    # convert to °C
    Tin_range = Tin_range - 273.15
    T_range = T_range - 273.15

    # generate a multiplicity plot
    plt.figure(1) 
    plt.plot(Tin_range, T_range,'k-')
    plt.axvline(T_in - 273.15, c='b')
    plt.xlabel("Inlet Temperature (°C)")
    plt.ylabel("Outlet Temperature (°C)")
    plt.savefig('example_6_6_4_multiplicity_plot.png')
    plt.savefig('example_6_6_4_multiplicity_plot.pdf')
    plt.show()
    
    return

# calculate the deliverables
if __name__=="__main__":
    deliverables()
