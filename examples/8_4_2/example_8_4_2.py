"""Calculations for Example 8.4.2 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ivodes

# global constants available to all functions
# given
k0_1 = 1.83E12 # L/mol/h
E_1 = 18000.0 # cal/mol
k0_2 = 5.08E13 # L/mol/h
E_2 = 20500.0 # cal/mol
dH_1 = -9000.0 # cal/mol
dH_2 = -7800.0 # cal/mol
Cp = 863.0 # cal/L/K
CA_0 = 2.0 # mol/L
CB_in = 0.5 # mol/L
T_0 = 40 + 273.15 # K
T_in = T_0
P = 1.0 # atm
V_A = 2000.0 # L
V_B = 8000.0 # L
t_f = 8.0 # h
t_1_range = np.array([1, 3, 5, 7])
# known
Re = 1.987 # cal/mol/K
Rw = 0.08206 # L-atm/mol/K
# calculated
nA_0 = CA_0*V_A

# define a global variables for the current value of t_1 and the current protocol stage
g_t_1 = float('nan')
g_stage = -1

# SBSTR reactor function
def sbstr_model_variables(t_1):
    # make t_1 available to the derivatives function
    global g_t_1
    g_t_1 = t_1

    # set the current operational protocol stage
    global g_stage
    g_stage = 1

    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nA_0, 0, 0, 0, T_0, V_A])

    # set the stopping criterion
    f_var = 0
    f_val = t_1

    # solve the IVODE design equations
    t1, dep1, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                ,sbstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f'SBSTR model stage 1 issue: {message}')
        print('')
        input('Press return to continue.')
    
    # repeat for the second stage
    g_stage = 2
    ind_0 = t1[-1]
    dep_0 = np.array([dep1[0,-1], dep1[1,-1], dep1[2,-1]
                , dep1[3,-1], dep1[4,-1], dep1[5,-1]])
    f_val = t_f
    t2, dep2, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                ,sbstr_derivatives, odes_are_stiff=False)    
    if not success:
        print('')
        print(f'SBSTR model stage 2 issue: {message}')
        print('')
        input('Press return to continue.')
    
    # concatenate the stages
    t = np.concatenate([t1, t2])
    nA = np.concatenate((dep1[0,:], dep2[0,:]))
    nB = np.concatenate((dep1[1,:], dep2[1,:]))
    nD = np.concatenate((dep1[2,:], dep2[2,:]))
    nU = np.concatenate((dep1[3,:], dep2[3,:]))
    T = np.concatenate((dep1[4,:], dep2[4,:]))
    V = np.concatenate((dep1[5,:], dep2[5,:]))

    # return the reactor model variables
    return t, nA, nB, nD, nU, T, V

# SBSTR derivatives function
def sbstr_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nD = dep[2]
    nU = dep[3]
    T = dep[4]
    V = dep[5]

    # calculate the additional unknowns
    if g_stage == 1:
        Vdot_in = V_B/g_t_1
    else:
        Vdot_in = 0.0
    k_1 = k0_1*np.exp(-E_1/Re/T)
    k_2 = k0_2*np.exp(-E_2/Re/T)
    CA = nA/V
    CB = nB/V
    r_1 = k_1*CA*CB
    r_2 = k_2*CB**2
    nDotB_in = Vdot_in*CB_in

    # evaluate the derivatives
    dnAdt = -V*r_1
    dnBdt = nDotB_in - V*(r_1 + 2*r_2)
    dnDdt = V*r_1
    dnUdt = V*r_2
    dTdt = (-Vdot_in*Cp*(T - T_in) - V*r_1*dH_1 - V*r_2*dH_2 + P*Vdot_in*Re/Rw)/(V*Cp)
    dVdt = Vdot_in

    # return the derivatives
    return [dnAdt, dnBdt, dnDdt, dnUdt, dTdt, dVdt]

# deliverables function
def deliverables():
    # allocate storage for the results
    f_B = np.ones_like(t_1_range)*float('nan')
    S_DU = np.ones_like(t_1_range)*float('nan')
    Y_DB = np.ones_like(t_1_range)*float('nan')

    # loop through the t_1 values
    for i, t_1 in enumerate(t_1_range):
        # solve the design equations
        t, nA, nB, nD, nU, T, V = sbstr_model_variables(t_1)

        # calculate and save the quantities of interest
        f_B[i] = 100*(V_B*CB_in - nB[-1])/(V_B*CB_in)
        S_DU[i] = nD[-1]/nU[-1]
        Y_DB[i] = 100*nD[-1]/((V_B*CB_in))
    
    # tabulate, show and save the results
    results_df = pd.DataFrame({'Stage 1 Duration (h)': t_1_range, 'Conversion (%)':f_B
                , 'Selectivity D/U': S_DU, 'Yield D-B (%)': Y_DB})
    print('')
    print(results_df)
    print('')
    results_df.to_csv('example_8_4_2_results.csv',index=False)

# execution command
if __name__ == '__main__':
    deliverables()
    