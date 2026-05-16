"""Calculations for the Class 18 Learning Activity from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
CA_0 = 2 # mol /L
V = 0.4 # L
Cp = 0.42 # cal /g /°C
rho = 879 # g /L
k01 = 2.5E8 # L /mol /min
E1 = 14300 # cal /mol
dH1 = -16600 # cal /mol
T0 = 20 + 273.15 # K
tf = 60 # min
# known
R = 1.987 # cal / mol / K
# calculated
nA0 = CA_0*V

# BSTR reactor function
def bstr_model_variables():
    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nA0, 0, T0])

    # set the stopping criterion
    f_var = 0
    f_val = tf

    # solve the design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , bstr_derivatives, odes_are_stiff=True)
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model issue: {message}")
        print('')
        input("Press Enter to continue.")
    
    # return the bstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:]

# BSTR derivatives function
def bstr_derivatives(t, dep):
    # extract the dependent variables
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]

    # calculate the additional unknowns
    k1 = k01 * np.exp(-E1/(R*T))
    CA = nA / V
    r1 = k1 * CA

    # evaluate the derivatives
    dnA_dt = -r1 * V
    dnZ_dt = r1 * V
    dT_dt = ( - r1*V*dH1) / (rho*V*Cp)

    # return the derivatives
    return np.array([dnA_dt, dnZ_dt, dT_dt])

# deliverables function
def deliverables():
    # solve the BSTR design equations
    t, nA, nZ, T = bstr_model_variables()

    # tabulate, show, and save results
    data = [["Conversion", 100*(nA0 - nA[-1])/nA0, "%"]
            ,["Maximum T", np.max(T) - 273.15, "°C"]]
    results_df = pd.DataFrame(data)
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_18_adiabatic_results.csv',index=False)

    plt.figure(1)
    plt.plot(t, T-273.15)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Temperature (°C)')
    plt.savefig('activity_18_adiabatic_T_profile.pdf')
    plt.show(block=False)

    plt.figure(2)
    plt.plot(t, 100*(nA0-nA)/nA0)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Conversion of A (%)')
    plt.ylim(bottom=0, top=100)
    plt.savefig('activity_18_adiabatic_fA_profile.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    