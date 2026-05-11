"""Calculations for the Class 15 Learning Activity in REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes
from reb_utils import solve_ates

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
CA_0 = 3.0/1000 # mol/cc
CB_0 = 3.0/1000 # mol/cc
T_0 = 50 + 273.15 # K
T_max = 90 + 273.15 # K
A = 66 # cm^2
V_ex = 40 # cc
U = 35*252.1*0.001076/60*1.8 # cal/cm^2/min/K
T_ex_in = 40 + 273.15 # K
t_turn = 30 # min
Cp = 0.35 # cal/g/K
rho = 0.93 # g/cc
k0_1 = 1.24E13*60 # cc/mol/min
E_1 = 20000 # cal/mol
dH_1 = -80000/4.184 # cal/mol
rho_ex = 1 # g/cc
Cp_ex = 1 # cal/g/K
VA_0 = 250 # cc
T_ex_0 = 40 + 273.15 # K
VB_0 = 250 # cc
fA = 0.9
# known
R = 1.987 # cal/mol/K
# calculated
V = VA_0 + VB_0
nA_0 = CA_0*VA_0
nB_0 = CB_0*VB_0
nA_f = nA_0*(1 - fA)

# define a global variable for the coolant flow rate
global g_m_ex
g_m_ex = float("NaN")

# BSTR reactor model function
def bstr_model_variables(mDot_ex):
    # make the coolant flow rate available to the derivatives function
    global g_m_ex
    g_m_ex = mDot_ex

    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nA_0, nB_0, 0, 0, T_0, T_ex_0])

    # define the stopping criterion
    f_var = 1
    f_val = nA_f

    # solve the design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var
            , f_val, bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model function issue: {message}")
        print('')
        input('Press return to continue.')
    
    # return the bstr reactor variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:]

# BSTR derivatives function
def bstr_derivatives(t, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    T = dep[4]
    Tex = dep[5]
    
    # calculate the additional unknowns
    k = k0_1*np.exp(-E_1/R/T)
    CA = nA/V
    CB = nB/V
    r = k*CA*CB
    Q = U*A*(Tex - T)

    # calc the derivatives
    dnAdt = -r*V
    dnBdt = -r*V
    dnYdt = r*V
    dnZdt = r*V
    dTdt = (Q - r*dH_1*V)/(rho*V*Cp)
    dTexdt = (-Q - g_m_ex*Cp_ex*(Tex - T_ex_in))/(rho_ex*V_ex*Cp_ex)

    # return the derivatives
    return np.array([dnAdt, dnBdt, dnYdt, dnZdt, dTdt, dTexdt])

# coupled unknown residual function
def coupled_unknown_residual(m_ex_guess):
    # solve the BSTR design equations using the guess
    t, nA, nB, nY, nZ, T, Tex = bstr_model_variables(m_ex_guess)

    # evaluate and return the residual
    epsilon = max(T) - T_max
    return epsilon

# deliverables function
def deliverables():
    # guess the coolant flow rate
    m_ex_guess = 100 # g/min

    # calculate the coolant flow rate
    soln, success, message = solve_ates(coupled_unknown_residual, m_ex_guess)
    mDot_ex = soln[0]

    # check for solver issues
    if not success:
        print('')
        print(f"Issue solving for the coupled unknown: {message}")
        print('')
        input('Press return to continue.')
    
    # solve the BSTR design equations
    t, nA, nB, nY, nZ, T, Tex = bstr_model_variables(mDot_ex)

    # calculate the net rate
    net_rate = nZ/(t + t_turn)

    # generate, show, and save the requested graph
    plt.figure(1)
    plt.plot(t, net_rate)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Net Rate (mol min$^{-1}$)')
    plt.title(f"{mDot_ex:.1f} g/min Coolant Flow")
    plt.savefig('activity_15_r_vs_t.png')
    plt.savefig('activity_15_r_vs_t.pdf')
    plt.show(block=False)

    # for discussion, solve using 90% and 110% of the base flow rate
    t90, nA, nB, nY, nZ, T90, Tex = bstr_model_variables(0.9*mDot_ex)
    net_rate_90 = nZ/(t90 + t_turn)
    t110, nA, nB, nY, nZ, T110, Tex = bstr_model_variables(1.1*mDot_ex)
    net_rate_110 = nZ/(t110 + t_turn)

    # plot, show and save net rate vs t
    plt.figure(2)
    plt.plot(t90, net_rate_90, label=f"{0.9*mDot_ex:.1f} g /min")
    plt.plot(t, net_rate, label=f"{mDot_ex:.1f} g /min")
    plt.plot(t110,net_rate_110, label=f"{1.19*mDot_ex:.1f} g /min")
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Net Rate (mol min$^{-1}$)')
    plt.legend()
    plt.title("Effect of Coolant Flow Rate (g/min)")
    plt.savefig('activity_15_r_vs_t_mex.png')
    plt.savefig('activity_15_r_vs_t_mex.pdf')
    plt.show(block=False)

    # plot, show and save T vs t
    plt.figure(3)
    plt.plot(t90, T90 - 273.15, label=f"{0.9*mDot_ex:.1f} g /min")
    plt.plot(t, T-273.15, label=f"{mDot_ex:.1f} g /min")
    plt.plot(t110,T110 - 273.15, label=f"{1.19*mDot_ex:.1f} g /min")
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Temperature (°C)')
    plt.legend()
    plt.title("Effect of Coolant Flow Rate (g/min)")
    plt.savefig('activity_15_T_vs_t.png')
    plt.savefig('activity_15_T_vs_t.pdf')
    plt.show()

    return

# execution command
if __name__ == '__main__':
    deliverables()
    