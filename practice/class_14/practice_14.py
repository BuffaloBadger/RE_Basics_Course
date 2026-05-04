"""Calculations for the Class 14 Practice Assignment from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
CB0 = 4.0 # mol/L
CA0 = 10.0 # mol/L
dH = -44000.0 # cal/mol
k0 = 8.11E12 # L/mol/s
E = 17700.0 # cal/mol
V = 25.0 # L
VA0 = 4.0 # L
T0 = 20 + 273.15 # K
Tex0 = 20 + 273.15 # K
TexIn = 20 + 273.15 # K
mEx = 1.0 # kg/min
Vex = 0.5 # L
VB0 = 10.0 # L
A = 0.6 # ft^2
U = 1.13E4 # cal/ft^2/h/K
rho = 1.0 # g/cm^3
Cp = 1.0 # cal/g/K
# convert to consistent units
CA0 = CA0/1000 # mol/cm^3
CB0 = CB0/1000 # mol/cm^3
k0 = k0*1000*60 # cm^3/mol/min
VA0 = VA0*1000 # cm^3
mEx = mEx*1000 # g/min
Vex = Vex*1000 # cm^3
VB0 = VB0*1000 # cm^3
V = VA0 + VB0;
UA = U*A/60 # cal/min/K
# known
R = 1.987
# calculated
nA0  = CA0*VA0
nB0 = CB0*VB0

# BSTR reactor function
def bstr_model_variables(tf):
    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nA0, nB0, 0, 0, T0, Tex0])

    # set the stopping criterion
    f_var = 0
    f_val = tf

    # solve the design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var
            ,f_val, bstr_derivatives, odes_are_stiff=True)
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR function issue: {message}")
        print('')
        input('Press return to continue.')
    
    # return the bstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:]

# BSTR derivatives function
def bstr_derivatives(ind,dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nS = dep[2]
    nW = dep[3]
    T = dep[4]
    Tex = dep[5]

    # calculate the additional unknowns
    k = k0*np.exp(-E/(R*T))
    CA = nA/V
    CB = nB/V
    r = k*CA*CB
    Q = UA*(Tex - T)

    # evaluate the derivatives
    dnAdt = -V*r
    dnBdt = -V*r
    dnSdt = V*r
    dnWdt = V*r
    dTdt = (Q - V*r*dH)/(rho*V*Cp)
    dTexdt = (-Q - mEx*Cp*(Tex-TexIn))/(rho*Vex*Cp)

    # return the derivatives
    return dnAdt, dnBdt, dnSdt, dnWdt, dTdt, dTexdt

# deliverables function
def deliverables():
    # set tf
    tf = 10 # min

    # solve the BSTR design equations
    t, nA, nB, nS, nW, T, Tex = bstr_model_variables(tf)

    # calculate the concentration of A
    CA = nA/V

    # generate, show, and save the requested graphs
    plt.figure(1)
    plt.plot(t, CA*1000)
    plt.xlabel('t (min)')
    plt.xlim(left=0)
    plt.ylabel('C$_A$ (M)')
    plt.ylim(bottom=0)
    plt.savefig('practice_14_CA_vs_t.png')
    plt.savefig('practice_14_CA_vs_t.pdf')
    plt.show(block=False)

    plt.figure(2)
    plt.plot(t, T - 273.15)
    plt.xlabel('t (min)')
    plt.xlim(left=0)
    plt.ylabel('T (°C)')
    plt.savefig('practice_14_T_vs_t.png')
    plt.savefig('practice_14_T_vs_t.pdf')
    plt.show(block=False)

    plt.figure(3)
    plt.plot(t, Tex - 273.15)
    plt.xlabel('t (min)')
    plt.xlim(left=0)
    plt.ylabel('T$_{ex}$ (°C)')
    plt.savefig('practice_14_Tex_vs_t.png')
    plt.savefig('practice_14_Tex_vs_t.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    