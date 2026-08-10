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

# allocate global storage
g_Vdot_ex = float('nan') # for the current coolant flow rate
g_stage = 0 # for the current operational stage

# BSTR reactor model
def bstr_model_variables(Vdot_ex):
    # set the current exchange fluid flow rate
    global g_Vdot_ex
    g_Vdot_ex = Vdot_ex

    # set the current operational stage to 1
    global g_stage
    g_stage = 1

	# set the initial values for the first operational stage
    ind_0 = 0.0
    dep_0 = np.array([nA_0, 0.0, T_0, Tex_0])

	# define the stopping criterion for the first operational stage
    f_var = 3
    f_val = T_1
    
	# solve the design equations for the first operational stage
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , bstr_derivatives, odes_are_stiff=False)
    
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

    # set the current operational stage to 2
    g_stage = 2

    # set the initial values for the second operational stage
    ind_0 = t[-1]
    dep_0 = np.array([nA[-1], nZ[-1], T[-1], Tex[-1]])

	# define the stopping criterion for the second operational stage
    f_var = 3
    f_val = T_f
    
	# solve the design equations for the second operational stage
    t2, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not(success):
        print('')
        print(f"Stage 2 BSTR model issue: {message}")
        print('')
        input('Press return to continue.')

    # extract the dependent variable profiles
    nA2 = dep[0,:]
    nZ2 = dep[1,:]
    T2 = dep[2,:]
    Tex2 = dep[3,:]

    # concatenate the bstr model variables for the two stages
    t = np.concatenate((t, t2))
    nA = np.concatenate((nA, nA2))
    nZ = np.concatenate((nZ, nZ2))
    T = np.concatenate((T, T2))
    Tex = np.concatenate((Tex, Tex2))

    # return the bstr model variables
    return t, nA, nZ, T, Tex

# derivatives function
def bstr_derivatives(ind, dep):
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
    if g_stage == 1:
        dTdt = (Qdot_shell + Qdot_coil - V*r_1*dH_1)/V/Cp
        dTexdt = -Qdot_shell/rho_ex/V_shell/Cp_ex
    else:
        dTdt = (Qdot_shell - V*r_1*dH_1)/V/Cp
        dTexdt = -(Qdot_shell + g_Vdot_ex*rho_ex*Cp_ex*(Tex-Tex_in))/rho_ex/V_shell/Cp_ex
    
	# return the derivatives
    return [dnAdt, dnZdt, dTdt, dTexdt]

# deliverables function
def deliverables():
    # choose a range of coolant flow rates
    #coolant_flows = np.linspace(100.0, 500.0, 100)
    coolant_flows = np.linspace(175.0, 200.0, 100)

    # allocate storage for the corresponding net rates
    net_rates = np.zeros(100)

    # calculate the net rate for each coolant flow rate
    for i in range(0,100):
        # solve the reactor design equations
        t, nA, nZ, T, T_ex = bstr_model_variables(coolant_flows[i])
        
        # calculate the net rate
        net_rates[i] = nZ[-1]/(t[-1] + t_turn)

    # find the coolant flow where the net rate is maximized
    i_max = np.argmax(net_rates)
    Vdot_max = coolant_flows[i_max]

    # solve the reactor design equations using the optimum coolant flow
    t, nA, nZ, T, T_ex = bstr_model_variables(Vdot_max)
    
    # calculate the conversion vs time at the optimum coolant flow rate
    pct_conversion = 100*(nA_0 - nA)/nA_0
    
    # tabulate, show, and save the results
    max_net_rate = nZ[-1]/(t[-1] + t_turn)
    data = [["Optimum Coolant Flow", f"{Vdot_max:.0f}", "cm^3 min^-1^"],
            ["Maximum Net Rate", f"{max_net_rate:.4f}", "mol min^-1^"]]
    results_df = pd.DataFrame(data, columns=['item','value','units'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('example_7_4_3_results.csv', index=False)

    # generate, show, and save the graphs
    plt.figure(1) # net rate vs coolant flow
    plt.plot(coolant_flows, net_rates)
    plt.xlabel("Coolant Flow (cm$^3$/min)")
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
