"""Calculations for the Class 18 Practice Assignment from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
V = 11.5 # L
CA0 = 3.6 # mol/L
T0 = 35 + 273.15 # K
A = 0.012 # m^2
U = 1485 # cal/m^2/min/K
Tex = 373.15 # K
dHvap = 540 # cal/g
k0f = 6.297E5 # L/mol/min
Ef = 12650 # cal/mol
k0r = 5.148E18 # L/mol/min
Er = 30950 # cal/mol
dH = -18300 # cal/mol
rho = 1000 # g/L
Cp = 1.0 # cal/g/K
t_f = 120. # min
# known
R = 1.987 # cal/mol/K
# calculated
nA0 = CA0 * V # mol

# BSTR reactor function
def bstr_model_variables():
    # set initial values
    ind_0 = 0
    dep_0 = np.array([nA0, 0, 0, T0])

    # set the stopping criterion
    f_var = 0
    f_val = t_f

    # solve the BSTR design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model function issue: {message}")
        print('')
        input('press return to continue.')
    
    # return the bstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:]

# BSTR derivatives function
def bstr_derivatives(ind, dep):
    # extract the dependent variables
    nA, nY, nZ, T = dep

    # calculate the additional unknowns
    kf = k0f*np.exp(-Ef/R/T)
    kr = k0r*np.exp(-Er/R/T)
    CA = nA/V
    CY = nY/V
    CZ = nZ/V
    r = kf*CA**2 - kr*CY*CZ
    Q = U*A*(Tex-T)

    # evaluate the derivatives
    dnAdt = -V*r
    dnYdt = V*r
    dnZdt = V*r
    dTdt = (Q - V*r*dH)/(V*rho*Cp)

    # return the derivatives
    return np.array([dnAdt, dnYdt, dnZdt, dTdt])

# deliverables function
def deliverables():
    # solve the BSTR design equations
    t, nA, nY, nZ, T = bstr_model_variables()

    # calculate the conversion and the condensate flow rate
    fA = 100*(nA0 - nA)/nA0
    m_cond = U*A*(Tex - T) / dHvap

    # plot the results
    plt.figure(1)
    plt.plot(t, fA)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Conversion')
    plt.ylim(bottom=0)
    plt.savefig('practice_18_fA_vs_t.pdf')
    plt.show(block=False)

    plt.figure(2)
    plt.plot(t, T - 273.15)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Temperature (°C)')
    plt.savefig('practice_18_T_vs_t.pdf')
    plt.show(block=False)

    plt.figure(3)
    plt.plot(t, m_cond)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Condensate Flow Rate (g/min)')
    plt.ylim(bottom=0)
    plt.savefig('practice_18_m_cond_vs_t.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    