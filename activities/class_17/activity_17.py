"""Calculations for The Class 17 Learning Activity from REB, The Course"""

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
V = 2.0 # m^3
yA0 = 0.5
yB0 = 0.5
T0 = 450 # K
P0 = 7 # atm
dH1_298 = -6870 # cal/mol
CpA = 7.5 # cal/mol/K
CpB = 8.5 # cal/mol/K
CpY = 12.1 # cal/mol/K
CpZ = 5.7 # cal/mol/K
k01 = 83 # m^3 /mol /h
E1 = 10200 # cal/mol
tTurn = 20/60 # h
# known
R = 1.987 # cal/mol
Rpv = 8.206E-5 # m^3 atm/mol/K
# calculated
nA0 = yA0*P0*V/Rpv/T0
nB0 = yB0*P0*V/Rpv/T0

# BSTR reactor function
def bstr_model_variables(t_f):
    # set initial values
    ind_0 = 0
    dep_0 = np.array([nA0, nB0, 0, 0, T0, P0])

    # set the stopping criterion
    f_var = 0
    f_val = t_f

    # solve the design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , bstr_derivatives, odes_are_stiff=False)
    
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
    k1 = k01*np.exp(-E1/R/T)
    CA = nA/V
    CB = nB/V
    r = k1*CA*CB
    dH1 = dH1_298 + (CpZ + CpY - CpA - CpB)*(T-298)

    # create the mass matrix
    massMatrix = np.zeros((6,6))

    # edit the rows corresponding to the mole balances
    massMatrix[0,0] = 1
    massMatrix[1,1] = 1
    massMatrix[2,2] = 1
    massMatrix[3,3] = 1

    # edit the row corresponding to the energy balance
    massMatrix[4,4] = nA*CpA + nB*CpB + nY*CpY + nZ*CpZ
    massMatrix[4,5] = -V*R/Rpv

    # edit the row corresponding to the ideal gas law
    massMatrix[5,0] = Rpv*T
    massMatrix[5,1] = Rpv*T
    massMatrix[5,2] = Rpv*T
    massMatrix[5,3] = Rpv*T
    massMatrix[5,4] = Rpv*(nA + nB + nY + nZ)
    massMatrix[5,5] = -V

    # create the right-hand side vector
    rhs = np.array([-r*V, -r*V, r*V, r*V, -r*dH1*V, 0])

    # calculate and return the derivatives
    return sp.linalg.solve(massMatrix, rhs)

# deliverables function
def deliverables():
    # choose a large reaction time
    tf = 2 # h

    # solve the BSTR design equations
    t, nA, nB, nY, nZ, T, P = bstr_model_variables(tf)

    # calculate corresponding conversion and net rate
    fA = 100*(nA0 - nA)/nA0
    rNet = nY/(t + tTurn)

    # find the maximum net rate
    iOpt = np.argmax(rNet)

    # calculate the optimum time and the net rate, conversion and temperature at that time
    tOpt = t[iOpt]
    rNetMax = rNet[iOpt]
    f_at_max = fA[iOpt]
    T_at_max = T[iOpt]

    # tabulate, show, and save the results
    data = [["Optimum Reaction time", f"{tOpt*60:.0f}", "min"]
            ,["Maximum Net Rate", f"{rNetMax:.1f}", "mol/h"]
            ,["Conversion", f"{f_at_max:.1f}", "%"]
            ,["Temperature", f"{T_at_max:.0f}", "K"]]
    results_df = pd.DataFrame(data,columns=("Item", "Value", "Units"))
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_17_results.csv',index=False)

    # for discussion, plot the net rate vs time
    plt.figure(1)
    plt.plot(t*60,rNet)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Net Rate (mol h$^{-1}$)')
    plt.ylim(bottom=0)
    plt.savefig('activity_17_rNet_vs_t.pdf')
    plt.show(block=False)

    # for discussion, plot the temperature vs time
    plt.figure(2)
    plt.plot(t*60,T)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Temperature (K)')
    plt.savefig('activity_17_T_vs_t.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    