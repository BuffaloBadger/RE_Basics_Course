"""Calculations for Example 9.6.1 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ivodes

# global constants available to all functions
# given
dH_1 = 44800 # J /mol
k0_1 = 7.22E6*60 # mol /amt^2 /cm^3 /min
E_1 = 84100 # J /mol
L = 10.0*12*2.54 # cm
D = 2.54 # cm
T_ex = 200 + 273.15 # K
U = 7.48E4/60/12**2/2.54**2 # J /min /cm^2 /K
yA_in = 0.6
yB_in = 0.4
Vdot_in = 282E3 # cm^3 /min
P = 2.5 # atm
T_in = 175 + 273.15 # K
CpA = 18.0*4.184 # J /mol /K
CpB = 12.25*4.184 # J /mol /K
CpZ = 21.2*4.184 # J /mol /K
# known
Re = 8.314 # J /mol /K
Rw = 82.06 # cm^3 atm /mol /K
# calculated
nA_in = yA_in*P*Vdot_in/Rw/T_in
nB_in = yB_in*P*Vdot_in/Rw/T_in

# PFR reactor function
def prf_model_variables():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_in, nB_in, 0.0, T_in])

	# define the stopping criterion
    f_var = 0
    f_val = L
     
	# solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , pfr_derivatives, odes_are_stiff=False)

    # check that a solution was found
    if not(success):
        print('')
        print(f"PFR reactor function issue: {message}")
        print('')
        input('Press return to continue or CTRL-C to exit')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nZ = dep[2,:]
    T = dep[3,:]

    # return the pfr model variables
    return z, nA, nB, nZ, T

# PFR derivatives function
def pfr_derivatives(ind, dep):
	# extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nZ = dep[2]
    T = dep[3]

	# calculate the additional unknowns
    k_1 = k0_1*np.exp(-E_1/Re/T)
    PA = nA/(nA + nB + nZ)*P
    PB = nB/(nA + nB + nZ)*P
    r_1 = k_1*PA*PB

	# evaluate the derivatives
    dnAdz = -np.pi*D**2/4*r_1
    dnBdz = -np.pi*D**2/4*r_1
    dnZdz = np.pi*D**2/4*r_1
    dTdz = (np.pi*D*U*(T_ex - T) - np.pi*D**2/4*r_1*dH_1)/(nA*CpA + nB*CpB + nZ*CpZ)
    
	# return the derivatives
    return dnAdz, dnBdz, dnZdz, dTdz

# deliverables function
def deliverables():
    # solve the reactor design equations
    z, nA, nB, nZ, T = prf_model_variables()

    # calculate the quantities of interest
    T_out = T[-1]
    f_B = 100*(nB_in - nB[-1])/nB_in

    # tabulate the results
    data = [['T out', T_out - 273.15, '°C'],
		 ['B Conversion', f_B, '%']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)
    print('')

    # save the results
    results_df.to_csv('example_9_6_1_results.csv', index=False)

# execution command
if __name__ == '__main__':
    deliverables()
    