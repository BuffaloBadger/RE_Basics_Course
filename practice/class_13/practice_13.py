"""Calculations for Practice Assignment 13 of REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ivodes
import matplotlib.pyplot as plt

# set the resolution for graphs
plt.rc('savefig', dpi =300)

# global constants available to all functions
# given
V = 17.7 # ft^3
CA_in = 0.125 # lbmol/ft^3
Vdot_in = 0.15 # ft^3/min
T_in = 300 + 458.67 # °R
Cp = 73.1 # BTU/ft^3/°F
k01 = 3.04e8 # ft^3/lbmol/min
E1 = 29100 # BTU/lbmol
dH1 = 32800 # BTU/lbmol
rho_ex = 81.2 # lb/ft^3
Cp_ex = 1.06 # BTU/lb/°F
mDotex = 5.0 # lb/min
Tex_in = 410 + 458.67 # °R
Tex_0 = 355 + + 458.67 # °R
Vex = 7.0 # ft^3
A = 21.5 # ft^2
U = 1.71 # BTU/ft^2/min/°F
fA_ss = 0.853
T_0 = 332 + 458.67 # °R
tf = 10*60 # min
# known
R = 1.986 # BTU/lbmol/°R
# calculated constants
Vdot = Vdot_in
nA_in = CA_in*Vdot_in
nA_0 = nA_in*(1-fA_ss)
nZ_0 = nA_in*fA_ss

# cstr model function
def cstr_model_variables():
    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nA_0, nZ_0, T_0, Tex_0])

    # define the stopping criterion
    f_var = 0
    f_val = tf

    # solve the cstr design equations
    #odes_are_stiff = False
    odes_are_stiff = True
    ind, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , cstr_derivatives, odes_are_stiff)
    
    # check for solver issues
    if not success:
        print('')
        print(f"CSTR model function issue: {message}")
        print('')
        input('Press return to continue.')
    # return the cstr model variables
    return ind, dep[0,:], dep[1,:], dep[2,:], dep[3,:] # t, nA, nZ, T, Tex

# cstr derivatives function
def cstr_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]
    Tex = dep[3]

    # calculate the additional unknowns
    k1 = k01*np.exp(-E1/(R*T));
    CA = nA/Vdot
    r1 = k1*CA**2
    Q = U*A*(Tex - T)

    # evaluate the derivatives
    dnAdt = Vdot/V*(nA_in - nA - r1*V)
    dnZdt = Vdot/V*(-nZ + r1*V)
    dTdt = 1/V/Cp*(-Vdot*Cp*(T - T_in) - r1*V*dH1 + Q)
    dTexdt = 1/rho_ex/Vex/Cp_ex*(-Q - mDotex*Cp_ex*(Tex - Tex_in));

    # return the derivatives
    return dnAdt, dnZdt, dTdt, dTexdt

# deliverables function
def deliverables():
    # get the cstr model variables
    t, nA, nZ, T, Tex = cstr_model_variables()

    # calculate the concentration of Z
    CZ = nZ/Vdot

    # generate, show, and save the graphs
    plt.figure(1)
    plt.plot(t,CZ)
    plt.xlabel('Time (min)')
    plt.ylabel('Concentration of Z (lbmol ft$^{-3}$)')
    plt.savefig('practice_13_CZ_vs_t.png')
    plt.savefig('practice_13_CZ_vs_t.pdf')
    plt.show(block=False)

    plt.figure(2)
    plt.plot(t,T - 458.67)
    plt.xlabel('Time (min)')
    plt.ylabel('Outlet Temperature (°F)')
    plt.savefig('practice_13_T_vs_t.png')
    plt.savefig('practice_13_T_vs_t.pdf')
    plt.show()

#execution command
if __name__=='__main__':
    deliverables()