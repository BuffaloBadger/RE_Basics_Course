"""Calculations for Example 8.4.3 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes, solve_ates

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
k0_1 = 1.192E15 # /min
E_1 = 97600 # J/mol
dH_1 =-58615 # J/mol
V_ex = 300 # cm^3
Tex_in = 60 + 273.15 # K
Vdot_ex = 250 # cm^3/min
UA = 260*4.184 # J/min/K
P = 1.0 # atm
VW_0 = 67 # cm^3
VZ_0 = 283 # cm^3
Vacid = 0.3 # cm^3
T_0 = 60 + 273.15 # K
Tex_0 = 60 + 273.15 # K
Cp = 2.68 # J/cm^3/K
V_A = 350 # cm^3
T_in = 21 + 273.15 # K
T_f = 65 + 273.15 # K
T_max = 95 + 273.15 #K
rho_A = 1.082 # g/cm^3
rho_W = 1.0 # g/cm^3
rho_Z = 1.049 # g/cm^3
M_A = 102 # g/mol
M_W = 18 # g/mol
M_Z = 60 # g/mol
Cp_A = 168.2 # J/mol/K
rho_ex = 1.0 # g/cm^3
Cp_ex = 1.0*4.184 # cal/g/K
# known
Re = 1.987*4.184 # J/mol/K
Rw = 82.06 # cm^3-atm/mol/K
# calculated
V_0 = VW_0 + VZ_0 + Vacid
m_ex = Vdot_ex*rho_ex
nW_0 = VW_0*rho_W/M_W
nZ_0 = VZ_0*rho_Z*M_Z

# global variables for the current inlet volumetric flow rate and protocol stage
g_Vdot_in = float('nan')

# SBSTR reactor function
def sbstr_model_variables(Vdot_in):
    # set the current value of the inlet volumetric flow rate
    global g_Vdot_in
    g_Vdot_in = Vdot_in

    # set the initial values for the first stage
    ind_0 = 0.0
    dep_0 = np.array([0.0, nW_0, nZ_0, T_0, Tex_0, V_0])

    # set the stopping criterion for the first stage
    t_1 = V_A/Vdot_in
    f_var = 0
    f_val = t_1

    # solve the design equations for the first stage
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , sbstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f'SBSTR model first stage issue: {message}')
        print('')
        input('Press return to continue.')
    
    # extract the dependent variables
    nA = dep[0,:]
    nW = dep[1,:]
    nZ = dep[2,:]
    T = dep[3,:]
    Tex = dep[4,:]
    V = dep[5,:]

    # repeat for the second stage
    g_Vdot_in = 0
    ind_0 = t[-1]
    dep_0 = np.array([dep[0,-1], dep[1,-1], dep[2,-1]
        , dep[3,-1], dep[4,-1], dep[5,-1]])
    f_var = 4
    f_val = T_f
    t2, dep2, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , sbstr_derivatives, odes_are_stiff=False)
    if not success:
        print('')
        print(f'SBSTR model second stage issue: {message}')
        print('')
        input('Press return to continue.')
    
    # combine the results for the two stages
    t = np.concatenate((t, t2))
    nA = np.concatenate((nA, dep2[0,:]))
    nW = np.concatenate((nW, dep2[1,:]))
    nZ = np.concatenate((nZ, dep2[2,:]))
    T = np.concatenate((T, dep2[3,:]))
    Tex = np.concatenate((Tex, dep2[4,:]))
    V = np.concatenate((V, dep2[5,:]))

    # return the sbstr model variables
    return t, nA, nW, nZ, T, Tex, V

# SBSTR derivatives function
def sbstr_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nW = dep[1]
    nZ = dep[2]
    T = dep[3]
    Tex = dep[4]
    V = dep[5]

    # calculate the additional unknowns
    k_1 = k0_1*np.exp(-E_1/Re/T)
    CA = nA/V
    r_1 = k_1*CA
    Qdot = UA*(Tex-T)
    nA_in = g_Vdot_in*rho_A/M_A
    
    # evaluate the derivatives
    dnAdt = nA_in - V*r_1
    dnWdt = - V*r_1
    dnZdt = V*r_1
    dTdt = (Qdot - nA_in*Cp_A*(T - T_in) - V*r_1*dH_1 
            + P*g_Vdot_in*Re/Rw)/(Cp*V)
    dTexdt = (-Qdot - m_ex*Cp_ex*(Tex-Tex_in))/(rho_ex*V_ex*Cp_ex)
    dVdt = g_Vdot_in

    # return the derivatives
    return [dnAdt, dnWdt, dnZdt, dTdt, dTexdt, dVdt]

# coupled unknown residual function
def coupled_unknown_residual(guess):
    # solve the SBSTR design equations using the guess
    t, nA, nW, nZ, T, Tex, V = sbstr_model_variables(guess[0])

    # return the residual
    return T_max - max(T)

# items of interest function
def items_of_interest():
    # define a guess for the feed rate
    feedGuess = 35

    # solve the implicit equation for Vdot
    soln, success, message = solve_ates(coupled_unknown_residual, feedGuess)

    # check for solver issues
    if not success:
            print('')
            print(f'ATE solver issue: {message}')
            print('')
            input('Press return to continue.')

    # extract the results
    Vdot_opt = soln[0]

    # solve the reactor design equations using the optimum feed rate
    t, nA, nW, nZ, T, Tex, V = sbstr_model_variables(Vdot_opt)

    # calculate the conversion
    f_A = 100*(V_A*rho_A/M_A - nA[-1])/(V_A*rho_A/M_A)

    # tabulate results
    data = [["Optimum Feed Rate",f"{Vdot_opt}","cm^3^ min^-1^"]
            ,["Processing Time",f"{t[-1]}","min"]
            ,["Acetic Anhydride Conversion",f"{f_A}","%"]]
    results_df = pd.DataFrame(data, columns=['item','value','units'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('example_8_4_3_results.csv',index=False)

    # create, display and save graphs
    plt.figure(1) 
    plt.plot(t, T-273.15)
    plt.xlabel("Time (min)")
    plt.xlim(left=0)
    plt.ylabel("Temperature (°C)")
    plt.savefig('example_8_4_3_T_vs_t.png')
    plt.show(block=False)

    # calculate the concentration profile
    CA = nA/V

    plt.figure(2) 
    plt.plot(t, CA)
    plt.xlabel("Time (min)")
    plt.xlim(left=0)
    plt.ylabel("Acetic Anhydride Concentration (mol cm$^{-3}$)")
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.savefig('example_8_4_3_CA_vs_t.png')
    plt.show()

    # as a check, calculate the processing time and maximum temperature for
    # volumetric flow rates slightly smaller and larger than the optimum

    Vdot = [0.95*Vdot_opt, 1.05*Vdot_opt]
    Tmax = [float('nan'), float('nan')]
    tProc = [float('nan'), float('nan')]

    for i, Vflow in enumerate(Vdot):
        t, nA, nW, nZ, T, Tex, V = sbstr_model_variables(Vflow)
        tProc[i] = t[-1]
        Tmax[i] = max(T) - 273.15

    # tabulate, show, and save the results
    check_df = pd.DataFrame({'Vdot':Vdot, 'Tmax':Tmax, 'tProc':tProc})
    print('')
    print(check_df)
    print('')
    check_df.to_csv('example_8_4_3_check.csv',index=False)

# execution command
if __name__ == '__main__':
    items_of_interest()
    