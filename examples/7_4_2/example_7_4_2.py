"""Calculations for Example 7.4.2 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates
from reb_utils import solve_ivodes

# global constants available to all functions
# given
V = 10.0E3 # cm^3
Vex = 1.4E3 # cm^3
Uex = 138. # cal /ft^2 /min /K
Aex = 1200./929. # ft^2
Tex_0 = 40. + 273.15 # K
Tex_in = 40. + 273.15 # K
mDot_ex = 100. # g /min
rho = 1.0 # g /cm^3
rho_ex = 1.0 # g /cm^3
Cp = 1.0 # cal /g /K
Cp_ex = 1.0 # cal /g /K
CA_0 = 5.0E-3 # mol /cm^3
CB_0 = 7.0E-3 # mol /cm^3
dH_1 = -16.7E3 # cal /mol
dH_2 = -14.3E3 # cal /mol
k0_1 = 9.74E12 # cm^3 /mol /min
E_1 = 20.1E3 # cal /mol
k0_2 = 2.38E13 # /min
E_2 = 25.3E3 # cal /mol
t_f = 30. # min
fA_f = 0.45
# known
R = 1.987 # cal /mol /K
# calculated
nA_0 = CA_0*V
nB_0 = CB_0*V
nA_f = nA_0*(1 - fA_f)

# BSTR model function
def bstr_model_variables(T_0):
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, nB_0, 0.0, 0.0, 0.0, T_0, Tex_0])

	# define the stopping criterion
    f_var = 0
    f_val = t_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                    , bstr_derivatives, odes_are_stiff=False)

    # check that a solution was found
    if not(success):
        print('')
        print(f"BSTR function issue: {message}")
        print('')
        input('Press return to continue.')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nX = dep[2,:]
    nY = dep[3,:]
    nZ = dep[4,:]
    T = dep[5,:]
    Tex = dep[6,:]

    # return the BSTR model variables
    return t, nA, nB, nX, nY, nZ, T, Tex

# BSTR derivatives function
def bstr_derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    T = dep[5]
    Tex = dep[6]

	# calculate the rate
    CA = nA/V
    CB = nB/V
    k_1 = k0_1*np.exp(-E_1/R/T)
    r_1 = k_1*CA*CB
    k_2 = k0_2*np.exp(-E_2/R/T)
    r_2 = k_2*CA

    # calculate the rate of heat exchange
    Qdot = Uex*Aex*(Tex - T)

	# evaluate the derivatives
    dnAdt = -(r_1 + r_2)*V
    dnBdt = -r_1*V
    dnXdt = r_1*V
    dnYdt = r_1*V
    dnZdt = r_2*V
    dTdt = (Qdot-(r_1*dH_1 + r_2*dH_2)*V)/rho/V/Cp
    dTexdt = (-Qdot - mDot_ex*Cp_ex*(Tex - Tex_in))/rho_ex/Vex/Cp_ex

	# return the derivatives
    return [dnAdt, dnBdt, dnXdt, dnYdt, dnZdt, dTdt, dTexdt]

# coupled unknown residual function
def coupled_unknown_residual(T0_guess):
    # solve the BSTR design equations
    t, nA, nB, nX, nY, nZ, T, Tex = bstr_model_variables(T0_guess)

    # extract the calculated molar amount
    nA_f_calc = nA[-1]

    # evaluate the residual
    resid = nA_f_calc - nA_f

    # return the residual
    return resid

# deliverables function
def deliverables():
    # set initial guess for T0
    T0_guess = Tex_in + 20

    # solve the implict equation for T0
    soln, success, message = solve_ates(coupled_unknown_residual, T0_guess)

    # check for solver issues
    if not(success):
        print('')
        print(f"Issue when solving for the coupled unknown: {message}")
        print('')
        input('Press return to continue.')

    # extract the result
    T_0 = soln[0]

    # solve the reactor design equations
    t, nA, nB, nX, nY, nZ, T, Tex = bstr_model_variables(T_0)

    # check the final values
    print('')
    print('Check:')
    print(f"final time: {t[-1]} min, should be 30")
    print(f"final conversion: {100*(nA_0 - nA[-1])/nA_0}%, should be 45%")

    # calculate the other quantities of interest
    T_0 = T_0 - 273.15
    T_f = T[-1] - 273.15
    Tex_out = Tex[-1] - 273.15
    sel_X_Z = nX[-1]/nZ[-1]

    # tabulate the results
    data = [['T_0', f'{T_0}', '°C'],
    ['Tf', f'{T_f}', '°C'],
    ['Te_out',f'{Tex_out}', '°C'],
    ['sel_X_Z',f'{sel_X_Z}','mol X per mol Z']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)
    print('')

    # save the results
    results_df.to_csv('example_7_4_2_results.csv', index=False)
    return

if __name__=="__main__":
    deliverables()
