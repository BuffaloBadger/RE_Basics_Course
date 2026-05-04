"""Calculations for Example 7.4.3 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for the graph
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
PA_0 = 1. # atm
PB_0 = 2. # atm
T_0 = 25 + 273.14 # K
V = 2. # L
Aex = 600. # cm^2
Uex = 0.6 # cal /cm^2 /min /K
Tex = 30 + 273.15 # K
k_0_1 =  3.34E9 # mol /cm^3 /min /atm^2
k_0_2 = 4.99E9 # mol /cm^3 /min /atm^2
E_1 = 20.5E3 # cal /mol
E_2 = 21.8E3 # cal /mol
dH_1 = -6300. # cal /mol
dH_2 = -6900. # cal /mol
Cp_A = 7.4 # cal /mol /K
Cp_B = 8.6 # cal /mol /K
Cp_D = 10.7 # cal /mol /K
Cp_Z = 5.2 # cal /mol /K
Cp_U = 10.3 # cal /mol /K
# known
Re = 1.987 # cal /mol /K
Rw = 82.057 # cm^3 atm /mol /K
# calculated
nA_0 = PA_0*V/Rw/T_0
nB_0 = PB_0*V/Rw/T_0
P_0 = PA_0 + PB_0

# BSTR model function
def bstr_model_variables(t_f):
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, nB_0, 0.0, 0.0, 0.0, T_0, P_0])

	# define the stopping criterion
    f_var = 0
    f_val = t_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                    , bstr_derivatives, odes_are_stiff=True)

    # check for solver issues
    if not(success):
        print('')
        print(f"BSTR model function issue: {message}")
        print('')
        input('Press return to continue')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nD = dep[2,:]
    nZ = dep[3,:]
    nU = dep[4,:]
    T = dep[5,:]
    P = dep[6,:]

    # return the reactor model variables
    return t, nA, nB, nD, nZ, nU, T, P

# derivatives function
def bstr_derivatives(ind, dep):
	# extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nD = dep[2]
    nZ = dep[3]
    nU = dep[4]
    T = dep[5]

	# calculate the rate
    PA = nA*Rw*T/V
    PB = nB*Rw*T/V
    PD = nD*Rw*T/V
    k_1 = k_0_1*np.exp(-E_1/Re/T)
    k_2 = k_0_2*np.exp(-E_2/Re/T)
    r_1 = k_1*PA*PB
    r_2 = k_2*PD*PB

    # calculate the rate of heat exchange
    Qdot = Uex*Aex*(Tex - T)

	# Create mass matrix, setting all elements to zero
    mass_matrix = np.zeros((7,7))

    # Add 1 on the diagonal for the first 5 rows
    mass_matrix[0,0] = 1.0
    mass_matrix[1,1] = 1.0
    mass_matrix[2,2] = 1.0
    mass_matrix[3,3] = 1.0
    mass_matrix[4,4] = 1.0

    # Add the elements for the energy balance
    mass_matrix[5,5] = nA*Cp_A + nB*Cp_B + nD*Cp_D + nZ*Cp_Z + nU*Cp_U
    mass_matrix[5,6] = -V*Re/Rw

    # Add the elements for the ideal gas law equation
    mass_matrix[6,0] = Rw*T
    mass_matrix[6,1] = Rw*T
    mass_matrix[6,2] = Rw*T
    mass_matrix[6,3] = Rw*T
    mass_matrix[6,4] = Rw*T
    mass_matrix[6,5] = Rw*(nA + nB + nD + nZ +nU)
    mass_matrix[6,6] = -V

    # Create right side vector
    rhs1 = -r_1*V
    rhs2 = (-r_1 -r_2)*V
    rhs3 = (r_1-r_2)*V
    rhs4 = (r_1+r_1)*V
    rhs5 = r_2*V
    rhs6 = Qdot -(r_1*dH_1 + r_2*dH_2)*V
    rhs7 = 0.0
    rhs = np.array([rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7])

    # Evaluate the derivatives
    derivs = sp.linalg.solve(mass_matrix, rhs)

    # Return the derivatives
    return derivs

# deliverables function
def deliverables():
    # choose a large final time
    tf = 60 # min

    # solve the BSTR design equations
    t, nA, nB, nD, nZ, nU, T, P = bstr_model_variables(tf)

    # calculate the corresponding yield
    yield_DA = nD/nA_0

    # find the reaction time where the yield is maximized
    i_opt = np.argmax(yield_DA)

    # find the corresponding reaction time
    t_opt = t[i_opt]

    # check that the optimum time is between 0 and t_f
    if (t_opt == 0) or (t_opt == tf):
        print('')
        print('The range of reaction times needs to be wider')
        print('')
        input('Press return to continue')

    # calculate the optimum yield and conversion
    yield_opt = yield_DA[i_opt]
    f_opt = 100*(nA_0 - nA[i_opt])/nA_0

    # tabulate the results
    data = [['t_opt', f'{t_opt}', 'min'],
    ['Y_D/A', f'{yield_opt}', 'mol D per initial mol A'],
    ['Conversion',f'{f_opt}', '%']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)
    print(' ')

    # save the results
    results_df.to_csv('example_7_4_3_results.csv' , index=False)
    
    # plot, show, and save the yield vs. the reaction time
    plt.figure(1)
    plt.plot(t, yield_DA)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Yield (mol D per initial mol A)")
    plt.xlim(left=0)
    plt.savefig('example_7_4_3_yield_vs_t.pdf')
    plt.savefig('example_7_4_3_yield_vs_t.png')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex3/example_7_4_3_yield_vs_t.png')
    plt.show()
    return

if __name__=="__main__":
    deliverables()
