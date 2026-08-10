"""Calculations for Example 9.6.5 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ivodes

# global constants available to all functions
# given
k01 = 6.632E8 # cm^6 /g /mol /s
E1 = 15000 # cal /mol
L = 1600 # cm
D = 2 # cm
rho_bed = 2.5 # g/cm^3
rho = 0.65 # g/cm^3
Cp = 0.5 # cal /g /K
Vdot_in = 0.228E3/60.0 # cm^3/s
CA_in = 7.84E-3 # mol/cm^3
CB_in = 2.32E-3 # mol/cm^3
T_in = 20 + 273.15 # K
dH1 = -17500 # cal/mol
# known
R = 1.987 # cal/mol/K
# calculated
nA_in = Vdot_in*CA_in
nB_in = Vdot_in*CB_in

# PFR reactor function
def pfr_model_variables():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_in, nB_in, 0.0, T_in])

	# define the stopping criterion
    f_var = 0
    f_val = L
     
	# solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , pfr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not(success):
        print('')
        print(f"PFR reactor model solver issue: {message}")
        print('')
        input('Press return to continue or CTRL-C to exit.')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nZ = dep[2,:]
    T = dep[3,:]

    # return the PFR model variables
    return z, nA, nB, nZ, T

# PFR derivatives function
def pfr_derivatives(ind, dep):
	# extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nZ = dep[2]
    T = dep[3]

	# calculate the additional unknowns
    Vdot = Vdot_in
    k_1 = k01*np.exp(-E1/R/T)
    CA = nA/Vdot
    CB = nB/Vdot
    r_1 = k_1*CA*CB

	# evaluate the derivatives
    dnAdz = -np.pi*D**2/4*rho_bed*r_1
    dnBdz = -np.pi*D**2/4*rho_bed*r_1
    dnZdz = np.pi*D**2/4*rho_bed*r_1
    dTdz = -np.pi*D**2/4*rho_bed*r_1*dH1/(Vdot_in*rho*Cp)

	# return the derivatives
    return dnAdz, dnBdz, dnZdz, dTdz

# deliverables function
def deliverables():
    # solve the PFR design equations
    z, nA, nB, nZ, T = pfr_model_variables()

    # calculate the quantities of interest
    fA = 100*(nA_in - nA[-1])/nA_in

    # tabulate, show, and save the results
    results = [['Conversion', f'{fA:.1f}', '%']
               ,['Temperature', f'{T[-1] - 273.15:.2f}', '°C']]
    results_df = pd.DataFrame(results, columns=['Item', 'Value', 'Units'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('example_9_6_5_results.csv', index=False)

# execution command
if __name__ == '__main__':
    deliverables()
    