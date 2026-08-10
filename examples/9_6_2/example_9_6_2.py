"""Calculations for Example 9.6.2 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates
from reb_utils import solve_ivodes

# global constants available to all functions
# given
Vpfr = 40E3 # cm^3
CA_in = 1E-3 # mol /cm^3
CB_in = 1.2E-3 # mol /cm^3
Vdot_in = 75E3 # cm^3 /min
k0_1 = 8.72E8 # cm^3 /mol /min
E_1 = 7200 # cal /mol
dH_1 = -10700 # cal /mol
Cp = 1.0 # cal /g /K
rho = 1.0 # g /cm^3
f_A = 0.95
# known
R = 1.987 # cal /mol /K
# calculated
nA_in = Vdot_in*CA_in
nB_in = Vdot_in*CB_in
nA_out = nA_in*(1-f_A)

# PFR reactor function
def pfr_model_variables(T_in):
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_in, nB_in, 0.0, 0.0, T_in])

	# define the stopping criterion
    f_var = 0
    f_val = Vpfr
     
	# solve the IVODEs
    V, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , pfr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not(success):
        print('')
        print(f"PFR model function issue: {message}")
        print('')
        input('Press return to continue or CTRL-C to exit')

    # extract the dependent variables
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]
    T = dep[4,:]

    # return the PFR model variables
    return V, nA, nB, nY, nZ, T

# PFR derivatives function
def pfr_derivatives(ind, dep):
	# extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nY = dep[2]
    nZ = dep[3]
    T = dep[4]

	# calculate the additional unknowns
    k_1 = k0_1*np.exp(-E_1/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r_1 = k_1*CA*CB

	# evaluate the derivatives
    dnAdV = -r_1
    dnBdV = -r_1
    dnYdV = r_1
    dnZdV = r_1
    dTdV = -r_1*dH_1/(Vdot_in*rho*Cp)

	# return the derivatives
    return dnAdV, dnBdV, dnYdV, dnZdV, dTdV

# coupled unknown residual function
def coupled_unknown_residual(guess):
    # extract the guess
    T_in = guess[0]

    # solve the reactor design equations
    V, nA, nB, nY, nZ, T = pfr_model_variables(T_in)

    # extract the calculated final value of nA
    nA_f = nA[-1]

    # evaluate the residual
    residual = nA_out - nA_f

    # return the residual
    return residual

# deliverables function
def deliverables():
    # set an initial guess for T_in
    initial_guess = 25 + 273.15

    # solve the implict equation for [missing constant]
    soln, success, message = solve_ates(coupled_unknown_residual,initial_guess)

    # check for solver issues
    if not(success):
        print('')
        print(f"Issue solving for the coupled unknown: {message}")
        print('')
        input('Press return to continue or CTRL-C to exit')

    # extract the result
    T_in = soln[0]

    # solve the reactor design equations
    V, nA, nB, nY, nZ, T = pfr_model_variables(T_in)

    # calculate the other quantities of interest
    T_f = T[-1]

    # tabulate the results
    data =[['Inlet T',T_in - 273.15,'°C'],['Outlet T',T_f - 273.15,'°C']]
    result_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(result_df)
    print('')

    # save the results
    result_df.to_csv('example_9_6_2_results.csv', index=False)

# execution command
if __name__ == '__main__':
    deliverables()
    