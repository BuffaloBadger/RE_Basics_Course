"""Calculations for the Class 23 Learning Activity from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
yA_in = 0.2
yB_in = 0.45
yC_in = 0.25
yI_in = 0.1
T_in = 220 + 273.15 # K
P_in = 3 # atm
Vdot_in = 1000 # cm^3 /min
D = 5 # cm
L = 1000 # cm
Tex = 185 + 273.15 # K
U = 1500 * 10E-4 # cal /cm^2 /min /K
CpA = 12.7 # cal /mol /K
CpB = 8.6
CpC = 11.3
CpX = 6.3
CpY = 14.4
CpZ = 10.8
CpI = 15.6
k0_1 = 1.0E5 # mol cm^-3 atm^-2
k0_2 = 5.0E5
E_1 = 19700 # cal /mol
E_2 = 21300
dH_1 = -28300
dH_2 = -29800
# known
Re = 1.987 # cal /mol /K
Rpv = 82.06 # cm^3 atm /mol /K
# calculated
nDotA_in = yA_in * P_in * Vdot_in / (Rpv*T_in)
nDotB_in = yB_in * P_in * Vdot_in / (Rpv*T_in)
nDotC_in = yC_in * P_in * Vdot_in / (Rpv*T_in)
nDotX_in = 0
nDotY_in = 0
nDotZ_in = 0
nDotI_in = yI_in * P_in * Vdot_in / (Rpv*T_in)
P = P_in

# PFR reactor function
def pfr_model_variables():
    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nDotA_in, nDotB_in, nDotC_in, nDotX_in, nDotY_in
        , nDotZ_in, nDotI_in, T_in])
    
    # set the stopping criterion
    f_var = 0
    f_val = L

    # solve the design equations
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , pfr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f'PFR model function issue: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit.')
    
    # return the pfr model variables
    return z, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:]\
        , dep[6,:], dep[7,:]

# PFR derivatives function
def pfr_derivatives(ind,dep):
    # extract the dependent variables
    nDotA = dep[0]
    nDotB = dep[1]
    nDotC = dep[2]
    nDotX = dep[3]
    nDotY = dep[4]
    nDotZ = dep[5]
    nDotI = dep[6]
    T = dep[7]

    # calculate the additional unknowns
    k_1 = k0_1*np.exp(-E_1/(Re*T))
    k_2 = k0_2*np.exp(-E_2/(Re*T))
    n_total = nDotA + nDotB + nDotC + nDotX + nDotY + nDotZ + nDotI
    PA = nDotA*P/n_total
    PB = nDotB*P/n_total
    PC = nDotC*P/n_total
    r_1 = k_1*PA*PB
    r_2 = k_2*PB*PC

    # evaluate the derivatives
    dnAdz = np.pi*D**2/4*(-r_1)
    dnBdz = np.pi*D**2/4*(-r_1 - r_2)
    dnCdz = np.pi*D**2/4*(-r_2)
    dnXdz = np.pi*D**2/4*(r_1 + r_2)
    dnYdz = np.pi*D**2/4*(r_1)
    dnZdz = np.pi*D**2/4*(r_2)
    dnIdz = 0
    denominator = nDotA*CpA + nDotB*CpB + nDotC*CpC + nDotX*CpX + nDotY*CpY\
        + nDotZ*CpZ + nDotI*CpI
    dTdz = (np.pi*D*U*(Tex - T) - np.pi*D**2/4*(r_1*dH_1 + r_2*dH_2))/denominator

    # return the design equation derivatives
    return [dnAdz, dnBdz, dnCdz, dnXdz, dnYdz, dnZdz, dnIdz, dTdz]

# deliverables function
def deliverables():
    # solve the PFR design equations
    z, nA, nB, nC, nX, nY, nZ, nI, T = pfr_model_variables()

    # calculate the quantities of interest
    fB = 100*(nDotB_in - nB[-1])/nDotB_in
    selectivity = nY[-1]/nZ[-1]
    T_out = T[-1] - 273.15

    # report the results
    print('')
    print(f'conversion: {fB:.2f}')
    print(f'selectivity: {selectivity:.2f}')
    print(f'T out: {T_out:.1f} °C')

    plt.figure(1)
    plt.plot(z,T-273.15)
    plt.xlabel('axial position, z')
    plt.ylabel('Temperature (°C)')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    