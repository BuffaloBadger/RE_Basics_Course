"""Calculations for the Class 17 Practice Assignment from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
k01 = 1.75E9 # L^0.5 mol^-1.5 s^-1
E1 = 35700 # cal/mol
dH1 = -12900 # cal/mol
CpA = 37.0 # cal/mol/K
CpB = 34.0 # cal/mol/K
CpY = 36.0 # cal/mol/K
CpZ = 38.0 # cal/mol/K
P0 = 1.0 # atm
T0 = 400 + 273.15 # K
yA0 = 0.5
yB0 = 0.5
fA = 0.5
V = 5.0 # L
# known
Rpv = 0.08206 # L-atm/mol/K
R = 1.987 # cal/mol/K
# calculated
n0 = P0*V/Rpv/T0
nA0 = yA0*n0
nB0 = yB0*n0
nAf = nA0*(1-fA)

# BSTR reactor function
def bstr_model_variables():
    # define the initial values
    ind_0 = 0
    dep_0 = np.array([nA0, nB0, 0, 0, T0, P0])

    # define the stopping criterion
    f_var = 1
    f_val = nAf

    # solve the BSTR design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , bstr_derivatives, odes_are_stiff=True)
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model function issue: {message}")
        print('')
        input('Press return to continue.')
    
    # return the bstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:]

# BSTR derivatives function
def bstr_derivatives(t, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nY = dep[2]
    nZ = dep[3]
    T = dep[4]
    P = dep[5]

    # calculate the additional unknowns
    k1 = k01*np.exp(-E1/(R*T))
    CA = nA/V
    CB = nB/V
    r1 = k1*CA*np.sqrt(CB)

    # form the mass matrix
    mm = np.zeros((6,6))

    # add the coefficients for the mole balances
    mm[0,0] = 1
    mm[1,1] = 1
    mm[2,2] = 1
    mm[3,3] = 1

    # add the coefficients for the energy balance
    mm[4,4] = nA*CpA + nB*CpB + nY*CpY + nZ*CpZ
    mm[4,5] = -V

    # add the coefficients for the differential gas law
    mm[5,0] = Rpv*T
    mm[5,1] = Rpv*T
    mm[5,2] = Rpv*T
    mm[5,3] = Rpv*T
    mm[5,4] = Rpv*(nA + nB + nY + nZ)
    mm[5,5] = -V

    # form the right side vector
    rhs = np.array([-r1*V, -r1*V, r1*V, r1*V, -r1*V*dH1, 0])

    # evaluate and return the derivatives
    derivs = sp.linalg.solve(mm, rhs)

    # Return the derivatives
    return derivs

# deliverables function
def deliverables():
    # solve the BSTR design equations
    t, nA, nB, nY, nZ, T, P = bstr_model_variables()

    # tabulate, show, and save the results
    data = [["Time", f"{t[-1]/60:.1f}", "min"]
            ,["Temperature", f"{T[-1] - 273.15:.0f}", "°C"]]
    results_df = pd.DataFrame(data,columns=["Item", "Value", "Units"])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('practice_17_results.csv',index=False)

    # check on pressure
    P_ratio_from_mole_balance = P[-1]/P0
    P_ratio_from_temperatures = T[-1]/T0

    print('')
    print(f"Pressure ratio from mole balance: {P_ratio_from_mole_balance:.2f}")
    print(f"Pressure ratio from temperatures: {P_ratio_from_temperatures:.2f}")

# execution command
if __name__ == '__main__':
    deliverables()
    