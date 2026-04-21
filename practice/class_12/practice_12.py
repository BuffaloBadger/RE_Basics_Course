"""Calculations for the Class 12 Practice Assignment of REB, The Course"""

# import libraries
import numpy as np
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
V = 1.5 # L
T_in = 50 + 273.15 # K
CA_in = 1.5 # mol/L
CB_in = 2.5 # mol/L
k01 = 1.4E14 # /min
E1 = 22000 # cal/mol
dH1 = -211000*0.239 # cal/mol
Cp = 1.3*1000 # cal/L/K
tau = 10 # min
T_0 = T_in
# known
R = 1.987 # cal/mol/K
# calculated
Vdot_in = V/tau
Vdot = Vdot_in
nA_in = CA_in*Vdot_in
nB_in = CB_in*Vdot_in

# cstr model function
def cstr_model_variables():
    # set the initial values
    ind_0 = 0.0
    dep_0 = np.array([0, 0, 0, 0, T_0])

    # set the stopping criterion
    f_var = 0
    f_val = 60 # min

    # solve the cstr design equations
    odes_are_stiff = False
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                            ,cstr_derivatives, odes_are_stiff)
    
    # check for solver issues
    if not success:
        print('')
        print('CSTR model function issue: ' + message)
        print('')
        input('Press return to continue.')

    # return the cstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:] # t, nA, nB, nY, nZ, T

# cstr derivatives function
def cstr_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nY = dep[2]
    nZ = dep[3]
    T = dep[4]

    # calculate the additional unknowns
    k1 = k01*np.exp(-E1/(R*T))
    CA = nA/Vdot
    r1 = k1*CA

    # evaluate the derivatives
    dnAdt = Vdot/V*(nA_in - nA - 2*r1*V)
    dnBdt = Vdot/V*(nB_in - nB - r1*V)
    dnYdt = Vdot/V*(-nY + 2*r1*V)
    dnZdt = Vdot/V*(-nZ + 2*r1*V)
    dTdt = -(Vdot*Cp*(T - T_in) + r1*V*dH1)/(V*Cp)

    # return the derivatives
    return np.array([dnAdt, dnBdt, dnYdt, dnZdt, dTdt])

# deliverables function
def deliverables():
    # solve the cstr design equations
    t, nA, nB, nY, nZ, T = cstr_model_variables()

    # generate, show, and save the outlet molar flow of Y and temperature vs. time
    plt.figure(1)
    plt.plot(t,nY)
    plt.xlabel('Time (min)')
    plt.ylabel('Outlet Molar Flow of Y (mol min$^{-1}$)')
    plt.savefig('practice_12_nY_vs_t.png')
    plt.savefig('practice_12_nY_vs_t.pdf')
    plt.show(block=False)

    plt.figure(2)
    plt.plot(t,T - 273.15)
    plt.xlabel('Time (min)')
    plt.ylabel('Temperature (°C)')
    plt.savefig('practice_12_T_vs_t.png')
    plt.savefig('practice_12_T_vs_t.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
