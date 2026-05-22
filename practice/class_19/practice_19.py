"""Calculations for the Class 19 Practice Assignment from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
CBin = 4.0E-3 # mol/cc
CA0 = 10.0E-3 # mol/cc
dH1 = -44000.0 # cal/mol
k01 = 8.11E15 # cc/mol/min
E1 = 17700.0 # cal/mol
V0 = 4.0E3 # cc
T0 = 20 + 273.15 # K
Tex0 = 20 + 273.15 # K
TexIn = 20 + 273.15 # K
mEx = 1.0E3 # g/min
Vex = 500 # cc
Aex = 0.6 # ft^2
U = 1.13E4/60 # cal/ft^2/min/K
rho = 1.0 # g/cm^3
Cp = 1.0 # cal/g/K
P = 1.0 # atm
Vdot_in = 50 # cc/min
VB0 = 10.0E3 # cc
Tin = 20 + 273.15 # K
# known
Ren = 1.987 # cal/mol/K
Rpv = 82.06 # cc-atm/mol/K
# calculated
nA0 = CA0*V0;

# SBSTR reactor function
def sbstr_model_variables():
    # define the initial values
    ind_0 = 0
    dep_0 = np.array([nA0, 0, 0, 0, T0, Tex0, V0])

    # define the stopping criterion
    f_var = 7
    f_val = V0 + VB0

    # solve the design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , sbstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f'SBSTR Model Function Issue: {message}')
        print('')
        input('Press return to continue or CTRL+C to exit')
    
    # return the sbstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:], dep[6,:]

# SBSTR derivatives function
def sbstr_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nC = dep[2]
    nD = dep[3]
    T = dep[4]
    Tex = dep[5]
    V = dep[6]

    # calculate the additional unknowns
    k1 = k01*np.exp(-E1/Ren/T)
    CA = nA / V
    CB = nB / V
    r1 = k1 * CA * CB
    nBin = CBin * Vdot_in
    Qdot = U*Aex*(Tex - T)

    # evaluate the derivatives
    dnA_dt = -V*r1
    dnB_dt = nBin - V*r1
    dnC_dt = V*r1
    dnD_dt = V*r1
    dT_dt = (Qdot - Vdot_in*Cp*rho*(T - Tin) - r1*V*dH1 + P*Vdot_in*Ren/Rpv)/(V*Cp*rho)
    dTex_dt = -(Qdot + mEx*Cp*(Tex - TexIn))/(rho*Vex*Cp)
    dV_dt = Vdot_in

    return np.array([dnA_dt, dnB_dt, dnC_dt, dnD_dt, dT_dt, dTex_dt, dV_dt])

# deliverables function
def deliverables():
    # solve the SBSTR design equations
    t, nA, nB, nC, nD, T, Tex, V = sbstr_model_variables()

    # calculate the concentration of A
    CA = nA / V

    # plot CA vs. t
    plt.figure(1)
    plt.plot(t, CA*1000)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Concentration of A (M)')
    plt.ylim(bottom=0)
    plt.savefig('practice_19_CA_vs_time.pdf')
    plt.show(block = False)

    # plot T vs. t
    plt.figure(2)
    plt.plot(t, T-273.15)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Temperature (°C)')
    plt.savefig('practice_19_T_vs_time.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    