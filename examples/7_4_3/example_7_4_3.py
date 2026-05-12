"""Calculations for Example 7.4.3 from REB, The Book"""

# import libraries
import numpy as np
from reb_utils import solve_ivodes
import pandas as pd
import matplotlib.pyplot as plt

# set plot resolution
plt.rc('savefig', dpi=300)

# global constants available in all functions
# given
k0_1 = 2.59E9 # /min
E_1 = 16500. # cal /mol
dH_1 = -22200. # cal /mol
CA_0 = 2. # mol /L
T_0 = 23  + 273.15 # K
Cp = 440. # cal /L /K
V = 4.0 # L
V_shell = 0.5 # L
A_shell = 0.6 # ft^2
U_shell = 1.13E4/60 # cal /ft^2 /min /K
Tex_in = 20 + 273.15 # K
rho_ex = 1.0 # g /cm^3
Cp_ex = 1.0 # cal /g /K
U_coil = 3.8E4/60 # cal /ft^2 /min /K
A_coil = 0.23 # ft^2
T_coil = 120 + 273.15 # K
Tex_0 = 23 + 273.15 # K
T_1 = 50 + 273.15 # K
T_f = 25 + 273.15 # K
t_turn = 25 # min
# known
R = 1.987 # cal /mol /K
# calculated
nA_0 = CA_0*V

# allocate global storage for the current coolant flow rate
g_Vdot_ex = float('nan')

# BSTR reactor model for the first stage of the protocal
def stage1_bstr_model_variables():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, 0.0, T_0, Tex_0])

	# define the stopping criterion
    f_var = 3
    f_val = T_1
    
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , stage1_bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not(success):
        print('')
        print(f"Stage 1 BSTR model issue: {message}")
        print('')
        input('Press return to continue.')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]
    Tex = dep[3,:]

    # return the first stage model variables
    return t, nA, nZ, T, Tex

# derivatives function for the first stage of operation
def stage1_bstr_derivatives(ind, dep):
	# extract the dependent variables
    nA = dep[0]
    T = dep[2]
    Tex = dep[3]

	# calculate the rate
    CA = nA/V
    k_1 = k0_1*np.exp(-E_1/R/T)
    r_1 = k_1*CA
    
    # calculate the rates of heat exchange
    Qdot_shell = U_shell*A_shell*(Tex - T)
    Qdot_coil = U_coil*A_coil*(T_coil - T)

	# evaluate the derivatives
    dnAdt = -V*r_1
    dnZdt = V*r_1
    dTdt = (Qdot_shell + Qdot_coil - V*r_1*dH_1)/V/Cp
    dTexdt = -Qdot_shell/rho_ex/V_shell/Cp_ex
    
	# return the derivatives
    return [dnAdt, dnZdt, dTdt, dTexdt]

# reactor model for the second stage of operation
def stage2_bstr_model_variables(t_0, nA_0, nZ_0, T_0, Tex_0):
	# set the initial values
    ind_0 = t_0
    dep_0 = np.array([nA_0, nZ_0, T_0, Tex_0])

	# define the stopping criterion
    f_var = 3
    f_val = T_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                , stage2_bstr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not(success):
        print('')
        print(f"Stage 2 BSTR model issue: {message}")
        print('')
        input('Press return to continue.')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]
    Tex = dep[3,:]

    # return the second stage model variables
    return t, nA, nZ, T, Tex

# derivatives function for the second stage of operation
def stage2_bstr_derivatives(ind, dep):
	# extract the dependent variables
    nA = dep[0]
    T = dep[2]
    Tex = dep[3]

	# calculate the rate
    CA = nA/V
    k_1 = k0_1*np.exp(-E_1/R/T)
    r_1 = k_1*CA

    # calculate the rates of heat exchange
    Qdot_shell = U_shell*A_shell*(Tex - T)

	# evaluate the derivatives
    dnAdt = -V*r_1
    dnZdt = V*r_1
    dTdt = (Qdot_shell - V*r_1*dH_1)/V/Cp
    dTexdt = -(Qdot_shell + g_Vdot_ex*rho_ex*Cp_ex*(Tex-Tex_in))/rho_ex/V_shell/Cp_ex
    
	# return the derivatives
    return [dnAdt, dnZdt, dTdt, dTexdt]

# deliverables function
def deliverables():
    # solve the reactor design equations for the first stage
    t_1, nA_1, nZ_1, T_1, Tex_1 = stage1_bstr_model_variables()

    # allow this function to set g_Vdot_ex
    global g_Vdot_ex

    # choose a range of coolant flow rates
    #coolant_flows = np.linspace(100.0, 500.0, 100)
    coolant_flows = np.linspace(175.0, 200.0, 100)

    # allocate storage for the corresponding net rates
    net_rates = np.zeros(100)

    # calculate the net rate for each coolant flow rate
    for i in range(0,100):
        # make the coolant flow rate globally available
        g_Vdot_ex = coolant_flows[i]

        # solve the reactor design equations
        t, nA, nZ, T, T_ex = stage2_bstr_model_variables(t_1[-1], nA_1[-1], nZ_1[-1]
                                                , T_1[-1], Tex_1[-1])
        
        # calculate the net rate
        net_rates[i] = nZ[-1]/(t[-1] + t_turn)

    # find the coolant flow where the net rate is maximized
    i_max = np.argmax(net_rates)
    Vdot_max = coolant_flows[i_max]

    # solve the reactor design equations using the optimum coolant flow
    g_Vdot_ex = Vdot_max
    t_2, nA_2, nZ_2, T_2, _ = stage2_bstr_model_variables(t_1[-1], nA_1[-1], nZ_1[-1]
                                                , T_1[-1], Tex_1[-1])
    
    # combine the profiles for the two stages
    t = np.concatenate((t_1, t_2))
    nA = np.concatenate((nA_1, nA_2))
    nZ = np.concatenate((nZ_1, nZ_2))
    T = np.concatenate((T_1, T_2))
    
    # calculate the conversion vs time at the optimum coolant flow rate
    pct_conversion = 100*(nA_0 - nA)/nA_0
    
    # tabulate, show, and save the results
    max_net_rate = nZ[-1]/(t[-1] + t_turn)
    data = [["Optimum Coolant Flow", f"{Vdot_max:.0f}", "g min^-1^"],
            ["Maximum Net Rate", f"{max_net_rate:.4f}", "mol min^-1^"]]
    results_df = pd.DataFrame(data, columns=['item','value','units'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('example_7_4_3_results.csv', index=False)

    # generate, show, and save the graphs
    plt.figure(1) # net rate vs coolant flow
    plt.plot(coolant_flows, net_rates)
    plt.xlabel("Coolant Flow (g/min)")
    plt.ylabel("Net Rate (mol/min)")
    plt.savefig('example_7_4_3_net_rate_vs_coolant_flow.png')
    plt.savefig('example_7_4_3_net_rate_vs_coolant_flow.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex3/example_7_4_3_net_rate_vs_coolant_flow.png')
    plt.show(block=False)

    plt.figure() # conversion profile
    plt.plot(t,pct_conversion)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Conversion (%)")
    plt.xlim(left=0)
    plt.ylim(bottom = 0)
    plt.savefig('example_7_4_3_conversion_profile.png')
    plt.savefig('example_7_4_3_conversion_profile.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex3/example_7_4_3_conversion_profile.png')
    plt.show(block=False)

    plt.figure() # temperature profile
    plt.plot(t,T - 273.15)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Temperature (°C)")
    plt.xlim(left=0)
    plt.savefig('example_7_4_3_temperature_profile.png')
    plt.savefig('example_7_4_3_temperature_profile.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex3/example_7_4_3_temperature_profile.png')
    plt.show()

    return

# execution command
if __name__=="__main__":
    deliverables()
