"""Calculations for Example 8.4.1 from REB, The Book"""

# import libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
CA_0 = 15.0 # mol/L
T_0 = 20 + 273.15 # K
CB_in = 15.0 # mol/L
T_in = T_0
V_0 = 5.0 # L
Vdot_in = 0.25/60 # L/s
V_B = 5.0 # L
V_ex = 2.5 # L
Tex_in = 20 + 273.15 # K
Vdot_ex = 2.5/60 # L/s
A = 2150 # cm^2
U = 73/60/60 # cal/cm^2/s/K
Tex_0 = 20 + 273.15 # K
rho = 1000.0 # g/L
Cp = 1.0 # cal/g/K
dH_1 = -13700.0 # cal/mol
k0_1 = 8.11E12 # L/mol/s
E_1 = 17700.0 # cal/mol
P = 1.0 # atm
# known
Re = 1.987 # cal/mol/K
Rw = 0.08206 # L-atm/mol/K
# calculated
nDotB_in = Vdot_in*CB_in
nA_0 = CA_0*V_0
mDot_ex = Vdot_ex*rho
V_f = V_0 + V_B

# SBSTR reactor function
def sbstr_model_variables():
    # set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, 0.0, 0.0, 0.0, T_0, Tex_0, V_0])

    # define the stopping criterion
    f_var = 7
    f_val = V_f
     
    # solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                , sbstr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not(success):
        print('')
        print(f"SBSTR function issue: {message}")
        print('')
        input('Press return to continue')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nS = dep[2,:]
    nW = dep[3,:]
    T = dep[4,:]
    Tex = dep[5,:]
    V = dep[6,:]

    # return all profiles
    return t, nA, nB, nS, nW, T, Tex, V

# SBSTR derivatives function
def sbstr_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nS = dep[2]
    nW = dep[3]
    T = dep[4]
    Tex = dep[5]
    V = dep[6]

    # calculate the additional unknowns
    k_1 = k0_1*np.exp(-E_1/Re/T)
    CA = nA/V
    CB = nB/V
    r = k_1*CA*CB
    Qdot = U*A*(Tex - T)

    # evaluate the derivatives
    dnAdt = -V*r
    dnBdt = nDotB_in -V*r
    dnSdt = V*r
    dnWdt = V*r
    dTdt = (Qdot - rho*Vdot_in*Cp*(T - T_in) - V*r*dH_1 + P*Vdot_in*Re/Rw) \
        /(rho*V*Cp)
    dTexdt = (-Qdot - mDot_ex*Cp*(Tex - Tex_in))/(rho*V_ex*Cp)
    dVdt = Vdot_in

    # return the derivatives
    return [dnAdt, dnBdt, dnSdt, dnWdt, dTdt, dTexdt, dVdt]

# deliverables function
def deliverables():
    # solve the reactor design equations
    t, nA, nB, nS, nW, T, Tex, V = sbstr_model_variables()

    # calculate the deliverables for plotting
    CA = nA/V
    T_C = T - 273.15
    t_min = t/60

    # tabulate, show and save the final concentration and temperature
    data =[["Final Concentration of A", CA[-1], "M"]
           ,["Final Temperature", T_C[-1], "°C"]]
    results_df = pd.DataFrame(data, columns=['Item', 'Value', 'Units'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('example_8_4_1_results.csv', index=False)
    results_df.to_csv('../../../RE_Basics/solutions/ch8_ex1/example_8_4_1_results.csv', index=False)
    
    # display and save the graphs
    plt.figure(1)
    plt.plot(t_min,CA)
    plt.xlabel("$Time \; (min)$")
    plt.xlim(left=0)
    plt.ylabel("$Concentration of A \; (mol \; L^{-1})$")
    plt.ylim(bottom=0)
    plt.savefig('example_8_4_1_CA_vs_t.png')
    plt.savefig('../../../RE_Basics/solutions/ch8_ex1/example_8_4_1_CA_vs_t.png')
    plt.savefig('example_8_4_1_CA_vs_t.pdf')
    plt.show(block=False)

    plt.figure(2)
    plt.plot(t_min,T_C)
    plt.xlabel("$Time \; (min)$")
    plt.xlim(left=0)
    plt.ylabel("$Temperature \; (°C)$")
    plt.savefig('example_8_4_1_T_vs_t.png')
    plt.savefig('../../../RE_Basics/solutions/ch8_ex1/example_8_4_1_T_vs_t.png')
    plt.savefig('example_8_4_1_T_vs_t.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    