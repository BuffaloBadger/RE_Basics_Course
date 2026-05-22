"""Calculations for Discussion of Example 8.4.2 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ivodes

# global constants available to all functions
# given
k0_1 = 1.83E12 # L/mol/h
E_1 = 18000.0 # cal/mol
k0_2 = 5.08E13 # L/mol/h
E_2 = 20500.0 # cal/mol
dH_1 = -9000.0 # cal/mol
dH_2 = -7800.0 # cal/mol
Cp = 863.0 # cal/L/K
CA_0 = 2.0 # mol/L
CB_in = 0.5 # mol/L
T_0 = 40 + 273.15 # K
T_in = T_0
P = 1.0 # atm
V_A = 2000.0 # L
V_B = 8000.0 # L
t_f = 8.0 # h
# known
Re = 1.987 # cal/mol/K
Rw = 0.08206 # L-atm/mol/K
# calculated
nA_0 = CA_0*V_A
nB_0 = V_B*CB_in
V = V_A + V_B

# BSTR model function
def bstr_model_variables():
    # set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, nB_0, 0.0, 0.0, T_0])

    # set the stopping criterion
    f_var = 0
    f_val = t_f
     
    # solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues   
    if not success:
        print('')
        print(f'BSTR model issue: {message}')
        print('')
        input('Press return to continue.')
    
    # return the bstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:]

# BSTR derivatives function
def bstr_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nD = dep[2]
    nU = dep[3]
    T = dep[4]

    # calculate the additional unknowns
    k_1 = k0_1*np.exp(-E_1/Re/T)
    k_2 = k0_2*np.exp(-E_2/Re/T)
    CA = nA/V
    CB = nB/V
    r_1 = k_1*CA*CB
    r_2 = k_2*CB**2

    # evaluate the derivatives
    dnAdt = -V*r_1
    dnBdt = -V*(r_1 + 2*r_2)
    dnDdt = V*r_1
    dnUdt = V*r_2
    dTdt = (-V*r_1*dH_1 - V*r_2*dH_2)/(V*Cp)

    # return the derivatives
    return [dnAdt, dnBdt, dnDdt, dnUdt, dTdt]

# deliverables function
def deliverables():
    # solve the BSTR design equations
    t, nA, nB, nD, nU, T = bstr_model_variables()

    # calculate the quantities of interest
    f_B = 100*(nB_0 - nB[-1])/nB_0
    S_DoverU = nD[-1]/nU[-1]
    Y_DfromB = 100*nD[-1]/nB_0
    
    # tabulate, display and save the results
    data = [["Conversion of B", f_B, "%"],
            ["Selectivity for D over U", S_DoverU," "],
            ["Yield of D from B", Y_DfromB, "%"]]
    results_df = pd.DataFrame(data, columns=["item","value","units"])
    print(results_df)
    results_df.to_csv("example_8_4_2_discussion.csv", index=False)

# execution command
if __name__ == '__main__':
    deliverables()
