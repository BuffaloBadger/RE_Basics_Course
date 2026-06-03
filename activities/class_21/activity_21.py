"""Calculations for the Class 21 Learning Activity from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
nDotTotal_in = 100 # mol/s
yIin = 0.2
yAin = 0.6
yBin = 0.2
Tin = 150 + 273.15 # K
k01 = 2E15 # L^0.5/mol-0.5/s
E1 = 25100 # cal/mol
CpA = 5 # cal/mol/K
CpB = 7 # cal/mol/K
CpY = 6.5 # cal/mol/K
CpZ = 5.7 # cal/mol/K
CpI = 4.2 # cal/mol/K
fB = 0.95;
dH1 = -10000 # cal/mol
P = 2 # atm
# known
Ren = 1.987 # cal /mol /K
Rpv = 82.057E-3 # L atm /mol K
# calculated
nDotAin = yAin*nDotTotal_in
nDotBin = yBin*nDotTotal_in
nDotIin = yIin*nDotTotal_in
nDotBout = nDotBin*(1-fB)

# PFR reactor function
def pfr_model_variables():
    # define initial values
    ind_0 = 0
    dep_0 = np.array([nDotAin, nDotBin, 0, 0, nDotIin, Tin])

    # define stopping criterion
    f_var = 2
    f_val = nDotBout

    # solve the design equations
    V, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            ,pfr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f'PFR model issue: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit')
    
    # return the pfr model variables
    return V, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:]

# PFR derivatives function
def pfr_derivatives(V, dep):
    # extract the dependent variables
    nDotA = dep[0]
    nDotB = dep[1]
    nDotY = dep[2]
    nDotZ = dep[3]
    nDotI = dep[4]
    T = dep[5]

    # calculate the additional unknowns
    k1 = k01 * np.exp(-E1 / (Ren * T))
    nDotTotal = nDotA + nDotB + nDotY + nDotZ + nDotI
    Vdot = nDotTotal*Rpv*T/P
    CA = nDotA/Vdot
    CB = nDotB/Vdot
    r1 = k1*CA*np.sqrt(CB)

    # evaluate the derivatives
    dnDotAdV = -r1
    dnDotBdV = -r1
    dnDotYdV = r1
    dnDotZdV = r1
    dnDotIdV = 0
    dTdV = -(dH1*r1)/(nDotA*CpA + nDotB*CpB + nDotY*CpY + nDotZ*CpZ + nDotI*CpI)

    return np.array([dnDotAdV, dnDotBdV, dnDotYdV, dnDotZdV, dnDotIdV, dTdV])

# deliverables function
def deliverables():
    # solve the PFR design equations
    V, nDotA, nDotB, nDotY, nDotZ, nDotI, T = pfr_model_variables()

    # tabulate, show and save the results
    results = [['Vpfr', f"{V[-1]:.2f}", 'L']
               ,['Outlet Temperature', f"{T[-1]-273.15:.2f}", '°C']]
    results_df = pd.DataFrame(results, columns=['Item', 'Value', 'Units'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_21_results.csv', index=False)

    # for discussion, plot the temperature profile
    plt.figure(1)
    plt.plot(V, T-273.15)
    plt.xlabel('Volume (L)')
    plt.xlim(left=0)
    plt.ylabel('Temperature (°C)')
    plt.savefig('activity_21_temperature_profile.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    