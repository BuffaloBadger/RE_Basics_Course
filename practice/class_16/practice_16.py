"""Calculations for the Class 16 Practice Assignment from REB, The Course"""

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
Ve = 1.4 # L
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
Tf = 35 + 273.15 # K
t_turn = np.array([0.25, 0.5, 0.75]) # h
# known
R = 1.987 # cal / mol / K
# calculated
nA0 = CA0 * V

# global variable indication the stage in the protocol
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
    # allocate storage for the maximum net rate and corresponding conversion
    rNet_max = np.ones_like(t_turn) * float('nan')
    conv_opt = np.ones_like(t_turn)*float('nan')

    # set up the graph
    plt.figure(1)
    plt.xlabel('Conversion (%)')
    plt.xlim(0, 100)
    plt.ylabel('Net Rate (moles Z h$^{-1}$)')

    # loop through the different values of t_turn
    for i, turnaround in enumerate(t_turn):
        # set a range of heating times
        t1_range = np.linspace(0.05, 1.5, 100)

        # allocate storate for the correponding plot data
        conv = np.ones_like(t1_range) * float('nan')
        rNet_Z = np.ones_like(t1_range) * float('nan')

        # loop through the different values of t1
        for j, t1 in enumerate(t1_range):
            # solve the BSTR design equations
            t, nA, nZ, T = bstr_model_variables(t1)

            # calculate the net rate of production of Z
            rNet_Z[j] = nZ[-1] / (t[-1] + turnaround)
        
            # calculate the conversion
            conv[j] = 100 * (nA0 - nA[-1]) / nA0

        # save the maximum net rate and corresponding conversion
        iMax = np.argmax(rNet_Z)
        conv_opt[i] = conv[iMax]
        rNet_max[i] = rNet_Z[iMax]

        # add the plot data for this value of t_turn to the graph
        plt.plot(conv, rNet_Z, label = f'{turnaround} h')

    # complete the graph
    plt.legend(title = 'Turnaround Time')
    plt.savefig('practice_16_rNet_vs_f.pdf')
    plt.show()

    # tabulate the maximum net rates
    results_df = pd.DataFrame({'turnaround time': t_turn, 'maximum net rate': rNet_max
            , 'corresponding conversion': conv_opt})
    print('')
    print(results_df)
    print('')
    results_df.to_csv('practice_16_results.csv',index=False)

# execution command
if __name__ == '__main__':
    deliverables()
    