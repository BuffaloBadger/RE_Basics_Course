"""Calculations for the Class 20 Practice Assignment from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates
from reb_utils import solve_ivodes

# global constants available to all functions
# given
T_0 = 65 + 273.15 # K
T_in = 65 + 273.15 # K
CA_0 = 10.0E-3 # mol /cc
CB_0 = 14.0E-3 # mol /cc
VA = 5.0E3 # cc
VB = 5.0E3 # cc
P = 1 # atm
Vdot_in_range = np.array([0.25, 0.5]) * 1.0E3 # cc /min
fA_f = 0.45
Vex = 1.4E3 # cm^3
U = 138. # cal /ft^2 /min /K
Aex = 1200./929. # ft^2
Tex_in = 40. + 273.15 # K
mDot_ex = 100. # g /min
rho = 1.0 # g /cm^3
rho_ex = 1.0 # g /cm^3
Cp = 1.0 # cal /g /K
Cp_ex = 1.0 # cal /g /K
Tex_0 = 40. + 273.15 # K
dH_1 = -16.7E3 # cal /mol
dH_2 = -14.3E3 # cal /mol
k0_1 = 9.74E12 # cm^3 /mol /min
E_1 = 20.1E3 # cal /mol
k0_2 = 2.38E13 # /min
E_2 = 25.3E3 # cal /mol
# known
R = 1.987 # cal /mol /K
Rpv = 82.057 # cm^3 atm /mol /K
# calculated
nA_0 = CA_0*VA
nB_0 = CB_0*VB
nA_f = nA_0 * (1 - fA_f)

# allocate global storage for the current volumetric feed rate
g_Vdot = float('nan')

# SBSTR model function
def sbstr_model_variables(Vdot):
    # set the feed rate for stage 1
    global g_Vdot
    g_Vdot = Vdot

	# set the initial values for stage 1
    ind_0 = 0.0
    dep_0 = np.array([0, nB_0, 0.0, 0.0, 0.0, T_0, Tex_0, VB])

	# define the stopping criterion for stage 1
    f_var = 8
    f_val = VA + VB
     
	# solve the IVODEs for stage 1
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                    , sbstr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not(success):
        print('')
        print(f"SBSTR Stage 1 issue: {message}")
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
    V = dep[7,:]

    # set the feed rate for stage 2
    g_Vdot = 0

    # set the initial values for stage 2
    ind_0 = t[-1]
    dep_0 = np.array([nA[-1], nB[-1], nX[-1], nY[-1], nZ[-1], T[-1], Tex[-1], V[-1]])

    # define the stopping criterion for stage 2
    f_var = 1
    f_val = nA_f

    # solve the IVODEs for stage 2
    t2, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                    , sbstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not(success):
        print('')
        print(f"SBSTR stage 2 issue: {message}")
        print('')
        input('Press return to continue.')
    
    # combine the results for the two stages
    t = np.concatenate([t, t2])
    nA = np.concatenate([nA, dep[0,:]])
    nB = np.concatenate([nB, dep[1,:]])
    nX = np.concatenate([nX, dep[2,:]])
    nY = np.concatenate([nY, dep[3,:]])
    nZ = np.concatenate([nZ, dep[4,:]])
    T = np.concatenate([T, dep[5,:]])
    Tex = np.concatenate([Tex, dep[6,:]])
    V = np.concatenate([V, dep[7,:]])

    # return the BSTR model variables
    return t, nA, nB, nX, nY, nZ, T, Tex, V

# BSTR derivatives function
def sbstr_derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    nX = dep[2]
    nY = dep[3]
    nZ = dep[4]
    T = dep[5]
    Tex = dep[6]
    V = dep[7]

	# calculate the additional unknowns
    nDot_A = g_Vdot*CA_0
    CA = nA/V
    CB = nB/V
    k_1 = k0_1*np.exp(-E_1/R/T)
    r_1 = k_1*CA*CB
    k_2 = k0_2*np.exp(-E_2/R/T)
    r_2 = k_2*CA
    Qdot = U*Aex*(Tex - T)

	# evaluate the derivatives
    dnAdt = nDot_A - (r_1 + r_2)*V
    dnBdt = -r_1*V
    dnXdt = r_1*V
    dnYdt = r_1*V
    dnZdt = r_2*V
    dTdt = (Qdot - rho_ex*g_Vdot*Cp_ex*(T - T_in) - (r_1*dH_1 + r_2*dH_2)*V
            + P*g_Vdot*R/Rpv)/rho/V/Cp
    dTexdt = (-Qdot - mDot_ex*Cp_ex*(Tex - Tex_in))/rho_ex/Vex/Cp_ex
    dVdt = g_Vdot

	# return the derivatives
    return [dnAdt, dnBdt, dnXdt, dnYdt, dnZdt, dTdt, dTexdt, dVdt]

# deliverables function
def deliverables():
    # allocate storage for the results
    tf = np.ones(2)*float('nan')
    Tf = np.ones_like(tf)*float('nan')
    Tex_f = np.ones_like(tf)*float('nan')
    sel_X_Z = np.ones_like(tf)*float('nan')

    # loop through the range of Vdot values
    for i, Vdot in enumerate(Vdot_in_range):
        # solve the reactor design equations
        t, nA, nB, nX, nY, nZ, T, Tex, V = sbstr_model_variables(Vdot)

        # calculate and save the quantities of interest
        tf[i] = t[-1]
        Tf[i] = T[-1] - 273.15
        Tex_f[i] = Tex[-1] - 273.15
        sel_X_Z[i] = nX[-1]/nZ[-1]

    # tabulate the results
    results_df = pd.DataFrame({'Feed Rate (L/min)': Vdot_in_range * 1.0E-3
            , 'Processing time (min)': tf
            ,'Final T (°C)': Tf, 'Final Coolant T (°C)': Tex_f
            , 'Selectivity (X/Z)': sel_X_Z})

    # display the results
    print(' ')
    print(results_df)
    print('')

    # save the results
    results_df.to_csv('practice_20_results.csv', index=False)
    return

if __name__=="__main__":
    deliverables()
