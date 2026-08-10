"""Calculations for the Class 20 Learning Activity from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
k01 = 2.193e+12 # gal/lbmol/h
E1 = 32400 # BTU/lbmol
k02 = 6.087e+14 # gal/lbmol/h
E2 = 36900 # BTU/lbmol
dH1 = -16190 # BTU/lbmol
dH2 = -14030 # BTU/lbmol
Cp = 7.2 # BTU/gal/°R
P = 1.0 # atm
V0 = 550. # gal
T0 = 100 + 459.7 # °R
CA0 = 0.015 # lbmol/gal
VB = 2000. # gal
Tin = 100 + 459.7 # °R
CBin = 0.004 # lbmol/gal
fB_min = 0.9
tTurn = 0.5 # h
# known
R = 1.987 # BTU/lbmol/°R
Rpv = 5.4584 # gal-atm/lbmol/°R
# calculated
nA0 = CA0*V0

# global variable for the current volumetric flow rate
g_Vdot = float('nan')

# SBSTR reactor function
def sbstr_model_variables(Vdot):
    # set the volumetric flow rate for stage 1
    global g_Vdot
    g_Vdot = Vdot

    # set the initial values for the first stage
    ind_0 = 0
    dep_0 = np.array([nA0, 0, 0, 0, T0, V0])

    # set the stopping criterion
    f_var = 6
    f_val = V0 + VB

    # solve the design equations for stage 1
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                , sbstr_derivatives, odes_are_stiff=False)
    
    # extract the dependent variables
    nA = dep[0,:]
    nB = dep[1,:]
    nD = dep[2,:]
    nU = dep[3,:]
    T = dep[4,:]
    V = dep[5,:]
    
    # check for solver issues
    if not success:
        print('')
        print(f'SBSTR Model stage 1 issue: {message}')
        print('')
        input('Press Enter to continue or CTRL+C to exit')
    
    # calculate the conversion of B during stage 1
    fB1 = 100*(VB*CBin - nB[-1])/VB*CBin

    # if the conversion is greater or equal to 90% there is no stage 2
    if fB1 >= 90:
        # return the model variables
        return t, nA, nB, nD, nU, T, V
    
    # set the volumetric flow rate for stage 2
    g_Vdot = 0

    # set the initial values equal to the final values from stage 1
    ind_0 = t[-1]
    dep_0 = np.array([nA[-1], nB[-1], nD[-1], nU[-1], T[-1], V[-1]])

    # set the stopping criterion for stage 2
    f_var = 2
    f_val = VB*CBin*(1 - fB_min)
    
    # solve the design equations for stage 2
    t2, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                , sbstr_derivatives, odes_are_stiff=False)
    
    # combine the results for the two stages
    t = np.concatenate([t, t2])
    nA = np.concatenate([nA, dep[0,:]])
    nB = np.concatenate([nB, dep[1,:]])
    nD = np.concatenate([nD, dep[2,:]])
    nU = np.concatenate([nU, dep[3,:]])
    T = np.concatenate([T, dep[4,:]])
    V = np.concatenate([V, dep[5,:]])
    
    # check for solver issues
    if not success:
        print('')
        print(f'SBSTR Model stage 2 issue: {message}')
        print('')
        input('Press Enter to continue or CTRL+C to exit')
    
    # return the reactor model variables
    return t, nA, nB, nD, nU, T, V

# SBSTR derivatives function
def sbstr_derivatives(ind, dep):
    # Extract the state variables
    nA = dep[0]
    nB = dep[1]
    nD = dep[2]
    nU = dep[3]
    T = dep[4]
    V = dep[5]

    # calculate the additional unknowns
    k1 = k01 * np.exp(-E1/(R*T))
    k2 = k02 * np.exp(-E2/(R*T))
    CA = nA / V
    CB = nB / V
    r1 = k1 * CA * CB
    r2 = k2 * CB**2
    nB_in = g_Vdot * CBin
    
    # Calculate the derivatives
    dnAdt = -V*r1
    dnBdt = nB_in - V*(r1 + 2*r2)
    dnDdt = V*r1
    dnUdt = V*r2
    dTdt = (-g_Vdot*Cp*(T-Tin)-V*r1*dH1 -V*r2*dH2 + P*g_Vdot*R/Rpv)/V/Cp
    dVdt = g_Vdot
    
    return np.array([dnAdt, dnBdt, dnDdt, dnUdt, dTdt, dVdt])

# deliverables function
def deliverables():
    # set a range of values for the volumetric flow rate
    #Vdot_range = np.linspace(20, 40, 100)
    Vdot_range = np.linspace(25, 30, 100)

    # allocate storage for the corresponding net rate and number of stages
    rNet = np.ones_like(Vdot_range)*float('nan')

    # for each value of the volumetric flow rate
    for i, Vdot in enumerate(Vdot_range):
        # solve the design equations
        t, nA, nB, nD, nU, T, V = sbstr_model_variables(Vdot)

        # calculate and save the corresponding net rate
        rNet[i] = nD[-1]/(t[-1] + tTurn)
    
    # find the index of the maximum net rate
    iMax = np.argmax(rNet)

    # tabulate, show, and save the results
    Vdot_opt = Vdot_range[iMax]
    max_rNet = rNet[iMax]
    data = [["Optimum Feed Rate", Vdot_opt, "gal/h"]
            ,["Net Rate", max_rNet, "lbmol/h"]]
    results_df = pd.DataFrame(data, columns=["Item", "Value", "Units"])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_20_results.csv',index=False)

    # plot the net rate vs. the feed rate
    plt.figure(1)
    plt.plot(Vdot_range, rNet)
    plt.xlabel('Feed Rate (gal/h)')
    plt.ylabel('Net Rate (lbmol/h)')
    plt.tight_layout()
    plt.savefig('activity_20_rNet_vs_VFR.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    