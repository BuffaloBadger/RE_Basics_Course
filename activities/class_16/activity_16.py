"""Calculations for the Class 16 Learning Activity from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
V = 10 # L
U1 = 76 # cal /h /cm^2 /K
A1 = 600 # cm^2
U2 = 91 # cal /h /cm^2 /K
A2 = 375 # cm^2
Tex_1 = 115 + 273.15 # K
Tex_2 = 30 + 273.15 # K
CA0 = 5 # mol / L
T0 = 25 + 273.15 # K
Cp = 1000 # cal /L /K
dH1 = 20000 # cal / mol
k01 = 2.85E15 # 1/h
E1 = 25000 # cal / mol
t1 = 1.0 # h
Tf = 35 + 273.15 # K
t_turn = 0.5 # h
# known
R = 1.987 # cal / mol / K
# calculated
nA0 = CA0 * V

# global variable indicating the stage in the protocol
g_stage = 1

# BSTR reactor function
def bstr_model_variables(t1):
    # set the protocol stage to 1
    global g_stage
    g_stage = 1

    # set the initial values for stage 1 of the protocal
    ind_0 = 0
    dep_0 = np.array([nA0, 0, T0])

    # set the stopping criterion for stage 1 of the protocol
    f_var = 0
    f_val = t1

    # solve the design equations for stage 1 of the protocol
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                ,bstr_derivatives, odes_are_stiff = False)
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR function stage 1 solver issue: {message}")
        print('')
        input("Press return to continue.")
    
    # set the protocol stage to 2
    g_stage = 2

    # set the initial values for stage 2 of the protocal
    ind_0 = t[-1]
    dep_0 = np.array([nA[-1], nZ[-1], T[-1]])

    # set the stopping criterion for stage 2 of the protocol
    f_var = 3
    f_val = Tf

    # solve the design equations for stage 2 of the protocol
    t2, dep2, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                ,bstr_derivatives, odes_are_stiff = False)
    nA2 = dep2[0,:]
    nZ2 = dep2[1,:]
    T2 = dep2[2,:]
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR function stage 2 solver issue: {message}")
        print('')
        input("Press return to continue.")

    # combine the results from stage 1 and stage 2
    t = np.concatenate((t, t2))
    nA = np.concatenate((nA, nA2))
    nZ = np.concatenate((nZ, nZ2))
    T = np.concatenate((T, T2))
    
    # return the bstr model variables
    return t, nA, nZ, T

# BSTR derivatives function
def bstr_derivatives(t, dep):
    # extract the dependent variables
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]

    # calculate the additional unknowns
    k1 = k01 * np.exp(-E1 / (R * T))
    CA = nA / V
    r1 = k1 * CA
    if g_stage == 1:
        U = U1
        A = A1
        Tex = Tex_1
    else:
        U = U2
        A = A2
        Tex = Tex_2
    Qdot = U * A * (Tex - T)

    # evaluate the derivatives
    dnA_dt = -r1 * V
    dnZ_dt = r1*V
    dT_dt = (Qdot - r1 * V * dH1) / (Cp * V)

    # return the derivatives
    return np.array([dnA_dt, dnZ_dt, dT_dt])

# deliverables function
def deliverables():
    # solve the BSTR design equations
    t, nA, nZ, T = bstr_model_variables(t1)

    # calculate the net rate of production of Z
    rNet_Z = nZ[-1] / (t[-1] + t_turn)

    # tabulate, show, and save the results
    data = [["Heating Time", f"{t1:.2f}", "h"]
            ,["Process Time", f"{t[-1]:.2f}", "h"]
            ,["Conversion", f"{100 * (nA0 - nA[-1]) / nA0:.2f}", "%"]
            ,["Net Rate", f"{rNet_Z:.2f}", "mol Z /h"]]
    results_df = pd.DataFrame(data, columns = ["Item", "Value", "Units"])
    print('')
    print(results_df)
    print('')
    results_df.to_csv("activity_16_results.csv", index = False)

    # plot the conversion and temperature profiles
    conv = 100 * (nA0 - nA) / nA0
    plt.figure(1)
    plt.plot(t, conv)
    plt.axvline(x = t1, color = 'red', linestyle = '--')
    plt.xlabel('Time (h)')
    plt.ylabel('Conversion (%)')
    plt.xlim(left = 0)
    plt.ylim(bottom = 0, top = 100)
    plt.savefig('activity_16_conversion_profile.pdf')
    plt.show(block = False)

    plt.figure(2)
    plt.plot(t, T - 273.15)
    plt.axvline(x = t1, color = 'red', linestyle = '--')
    plt.xlabel('Time (h)')
    plt.ylabel('Temperature (°C)')
    plt.xlim(left = 0)
    plt.savefig('activity_16_temperature_profile.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    