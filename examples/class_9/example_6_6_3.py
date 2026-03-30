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
# known
R = 8.314 # J /mol /K

# make T_in globally available
T_in = float('NaN')

# reactor model function
def reactor_model_variables(current_T_in, init_guess):
    # set T_in
    global T_in
    T_in = current_T_in
     
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
    # set a range for T_in
    T_in_range = np.linspace(75.,125.,100) + 273.15

    # allocate storage for conversion and temperature at each T_in value
    fA_range = np.zeros(100)
    T_range = np.zeros(100)

	# set the initial guess
    init_guess = np.array([nA_in, nB_in, 0.0, T_in_range[0] + 5.0])

    # calculate fA and T for each T_in
    for iT in range(0,100):

        # solve the reactor design equations
        solution = reactor_model_variables(T_in_range[iT], init_guess)

        # calculate the conversion
        fA_range[iT] = 100*(nA_in - solution[0])/nA_in
        T_range[iT] = solution[3]

        # use the result as the initial guess for the next T_in
        init_guess = solution

    # plot the results
    plt.figure(1) 
    plt.plot(T_in_range-273.15, fA_range)
    plt.xlabel("Inlet Temperature (°C)")
    plt.ylabel("Conversion of A (%)")

    # save and show the figure
    plt.savefig('example_6_6_3_fA_vs_Tin.png')
    plt.savefig('example_6_6_3_fA_vs_Tin.pdf')
    plt.show()

    # find the optimum T_in
    i_max = np.argmax(fA_range)

    # calculate the quantities of interest
    T_in_opt = T_in_range[i_max] - 273.15
    fA_max = fA_range[i_max]
    T_out = T_range[i_max] - 273.15

    # tabulate the results
    data = [['opt T_in',f'{T_in_opt}','°C'],['fA_max',f'{fA_max}','%']
         ,['T out',f'{T_out}','°C']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('example_6_6_3_results.csv',index=False)

    # check for a high temperature steady state at the lowest T_in
    init_guess = np.array([nA_in, nB_in, 0.0, T_in_range[0] + 400.0])
    solution = reactor_model_variables(T_in_range[0], init_guess)
    print(f"Outlet T at start of range with low guess: {T_range[0]} K")
    print(f"Outlet T at start of range with high guess: {solution[3]} K")

    return

if __name__=="__main__":
    deliverables()
