"""Calculations for Discussion of Example 9.6.3 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ivodes

# global constants available to all functions
# given
yA_in = 0.76
yI_in = 0.25
Vdot_in = 100 # cm^3 /s
mDot_in = 0.44 # g /s
P_in = 3.0 # atm
P = P_in
T_in = 400 + 273.15 # K
D = 2.5 # cm
L = 800. # cm
T_ex = 375 + 273.15 # K
U = 187E-1/3600 # J /s /cm^2 /K
Dp = 0.25 # cm
phi = 0.7
eps = 0.6
k0f = 9E17 # mol /cm^3 /s /atm
Ef = 285E3 # J /mol
k0r = 4.09E-4 # mol /cm^3 /s /atm^4
Er = 85E3 # J /mol
dH = 200E3 # J /mol
CpA = 11.7*4.184 # J /mol /K
CpY = 8.3*4.184 # J /mol /K
CpZ = 4.2*4.184 # J /mol /K
CpI = 5.8*4.184 # J /mol /K
mu = 0.027E-2 # g /cm /s
# known
Re = 8.314 # J /mol /K
Rw = 82.06 # cm^3 atm /mol /K
P_conv = 9.872E-7 # atm cm^2 /dyne
# calculated
G = mDot_in/(np.pi*D**2/4)
nA_in = yA_in*Vdot_in*P_in/Rw/T_in
nI_in = yI_in*Vdot_in*P_in/Rw/T_in

# PFR reactor function
def pfr_model_variables():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_in, 0.0, 0.0, nI_in, T_in])

	# define the stopping criterion
    f_var = 0
    f_val = L
     
	# solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
        , pfr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not(success):
        print('')
        print(f"PFR model function issue: {message}")
        print('')
        input('Press return to continue or CTRL-C to exit.')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nY = dep[1,:]
    nZ = dep[2,:]
    nI = dep[3,:]
    T = dep[4,:]

    # return the pfr model variables
    return z, nA, nY, nZ, nI, T

# PFR derivatives function
def pfr_derivatives(ind, dep):
	# extract the dependent variables
    nA = dep[0]
    nY = dep[1]
    nZ = dep[2]
    nI = dep[3]
    T = dep[4]

	# calculate the additional unknowns
    kf = k0f*np.exp(-Ef/Re/T)
    kr = k0r*np.exp(-Er/Re/T)
    ntot = nA + nY + nZ + nI
    PA = nA/ntot*P
    PY = nY/ntot*P
    PZ = nZ/ntot*P
    r = kf*PA - kr*PY*PZ**3

	# evaluate the derivatives
    dnAdz = -np.pi*D**2/4*r
    dnYdz = np.pi*D**2/4*r
    dnZdz = 3*np.pi*D**2/4*r
    dnIdz = 0.0
    dTdz = (np.pi*D*U*(T_ex - T) - np.pi*D**2/4*r*dH)/(nA*CpA + nY*CpY 
            + nZ*CpZ + nI*CpI)

	# return the derivatives
    return dnAdz, dnYdz, dnZdz, dnIdz, dTdz

# deliverables function
def deliverables():
    # solve the reactor design equations
    z, nA, nY, nZ, nI, T = pfr_model_variables()

    # calculate the quantities of interest
    T_out = T[-1] - 273.15
    fA = 100*(nA_in - nA[-1])/nA_in

    # tabulate the results
    data =[['Outlet T',T_out,'°C'], ['Conversion',fA,'%']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)
    print('')

    # save the results
    results_df.to_csv('example_9_6_3_discussion.csv',index=False)

# execution command
if __name__ == '__main__':
    deliverables()
    