"""Calculations for Example 9.6.5 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

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
T_in_range = np.array([20, 100]) + 273.15 # K
dH1 = -17500 # cal/mol
# known
R = 1.987 # cal/mol/K
# calculated
nA_in = Vdot_in*CA_in
nB_in = Vdot_in*CB_in

# PFR reactor function
def pfr_model_variables(T_in):
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
    # allocate storage for the conversion and temperature
    fA = np.ones_like(T_in_range)*float('nan')
    T_out = np.ones_like(T_in_range)*float('nan')

    # loop over the inlet temperature range
    for i, T_in in enumerate(T_in_range):
        # solve the PFR design equations
        z, nA, nB, nZ, T = pfr_model_variables(T_in)

        # calculate the quantities of interest
        fA[i] = 100*(nA_in - nA[-1])/nA_in
        T_out[i] = T[-1] - 273.15

    # plot the conversion and temperature
    plt.figure(1)
    plt.plot(T_in_range - 273.15, fA)
    plt.xlabel('Inlet Temperature (°C)')
    plt.ylabel('Conversion of A (%)')
    plt.ylim(bottom=0)
    plt.savefig('example_9_6_5_f_vs_Tin.pdf')
    plt.show(block=False)

    plt.figure(2)
    plt.plot(T_in_range - 273.15, T_out)
    plt.xlabel('Inlet Temperature (°C)')
    plt.ylabel('Outlet Temperature (°C)')
    plt.savefig('example_9_6_5_Tout_vs_Tin.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    