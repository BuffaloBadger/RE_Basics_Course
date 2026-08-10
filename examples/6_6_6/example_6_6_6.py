"""Calculations for Example 6.6.6 of REB, The Book"""

# import libraries
import numpy as np
from reb_utils import solve_ivodes
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi=300)

# constants available to all functions
# given
V = 500.0 # cm^3
Vdot_in = 1.0 # cm^3 /s
CA_in = 0.015 # mol /cm^3
CB_in = 0.015 # mol /cm^3
T_in = 50 + 273.15 # K
T_high = 264.75 + 273.15 # K
fA_high = 0.9749
T_low = 50.06 + 273.15 # K
fA_low = 0.03
T_u = 137.96 + 273.15 # K
fA_u = .3993
Cp = 0.35*4.184 # J /g /K
rho = 0.93 # g /cm^3
dH = -20000.0 # J /mol
k0 = 3.24E12 # cm^3 /mol /s
E = 105000.0 # J /mol
# known
R = 8.314 # J /mol /K
# calculated
nA_in = CA_in*Vdot_in
nB_in = CB_in*Vdot_in

# cstr model function
def cstr_model_variables(T_0, fA_0, t_f):
	# set the initial values
    ind_0 = 0.0
    extent = fA_0*nA_in
    dep_0 = np.array([nA_in - extent, nB_in - extent, extent, extent, T_0])

	# define the stopping criterion
    f_var = 0
    f_val = t_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , cstr_derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"  CSTR model function error: {message}")

    # extract the dependent variable cstr_model_variables
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]
    T = dep[4,:]

    # return all cstr_model_variables
    return t, nA, nB, nY, nZ, T

# cstr derivatives function
def cstr_derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    nY = dep[2]
    nZ = dep[3]
    T = dep[4]

	# calculate the rate
    k = k0*np.exp(-E/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r = k*CA*CB
    
	# evaluate the derivatives
    dnAdt = Vdot_in/V*(nA_in - nA - V*r)
    dnBdt = Vdot_in/V*(nB_in - nB - V*r)
    dnYdt = Vdot_in/V*( -nY + V*r)
    dnZdt = Vdot_in/V*( -nZ + V*r)
    dTdt = -(Cp*Vdot_in*rho*(T-T_in) + V*r*dH)/(V*Cp*rho)

	# return the derivatives
    return dnAdt, dnBdt, dnYdt, dnZdt, dTdt

# deliverables function
def deliverables():

    # set the duration of the transient
    t_f = 2500.0 # s

    # analyze a positive perturbation from the unsteady state
    dT = 1.0 # K
    [t_u_pos, nA, nB, nY, nZ, T_u_pos] = cstr_model_variables(T_u + dT, fA_u, t_f)

    # analyze a negative perturbation from the unsteady state
    [t_u_neg, nA, nB, nY, nZ, T_u_neg] = cstr_model_variables(T_u - dT, fA_u, t_f)

    # analyze a positive perturbation from the high temperature steady state
    dT = 20.0 # K
    [t_h_pos, nA, nB, nY, nZ, T_high_pos] = cstr_model_variables(T_high + dT, fA_high, t_f)

    # analyze a negative perturbation from the high temperature steady state
    [t_h_neg, nA, nB, nY, nZ, T_high_neg] = cstr_model_variables(T_high - dT, fA_high, t_f)

    # analyze a positive perturbation from the low temperature steady state
    [t_l_pos, nA, nB, nY, nZ, T_low_pos] = cstr_model_variables(T_low + dT, fA_low, t_f)

    # analyze a negative perturbation from the low temperature steady state
    [t_l_neg, nA, nB, nY, nZ, T_low_neg] = cstr_model_variables(T_low - dT, fA_low, t_f)

    # create, display and save the graphs
    plt.figure(1) # perturbations from the unsteady state
    plt.plot(t_u_pos/60.0,T_u_pos - 273.15, label = '+1 °C Perturbation Response')
    plt.plot(t_u_neg/60.0,T_u_neg - 273.15, label = '-1 °C Perturbation Response')
    plt.axhline(T_high - 273.15, c = 'k', ls = '-', label = 'Stable Steady State, S-high')
    plt.axhline(T_u - 273.15, c = 'k', ls = '--', label = 'Unstable Steady State, U-mid')
    plt.axhline(T_low - 273.15, c = 'k', ls = '-', label = 'Stable Steady State, S-low')
    plt.xlabel("Time (min)")
    plt.ylabel("Reacting Fluid Temperature (°C)")
    plt.savefig('example_6_6_6_response_at_U_mid.png')
    plt.savefig('example_6_6_6_response_at_U_mid.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch6_ex6/example_6_6_6_response_at_U_mid.png')
    plt.legend()
    plt.show()

    plt.figure(2) # perturbations from the high T state
    plt.plot(t_h_pos/60.0,T_high_pos - 273.15, label = '+20 °C Perturbation Response')
    plt.plot(t_h_neg/60.0,T_high_neg - 273.15, label = '-20 °C Perturbation Response')
    plt.axhline(T_high - 273.15, c = 'k', ls = '-', label = 'Stable Steady State, S-high')
    plt.xlabel("Time (min)")
    plt.ylabel("Reacting Fluid Temperature (°C)")
    plt.savefig('example_6_6_6_response_at_S_high.png')
    plt.savefig('example_6_6_6_response_at_S_high.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch6_ex6/example_6_6_6_response_at_S_high.png')
    plt.legend()
    plt.show()

    plt.figure(3) # perturbations from the low T state
    plt.plot(t_l_pos/60.0,T_low_pos - 273.15, label = '+20 °C Perturbation Response')
    plt.plot(t_l_neg/60.0,T_low_neg - 273.15, label = '-20 °C Perturbation Response')
    plt.axhline(T_low - 273.15, c = 'k', ls = '-', label = 'Stable Steady State, S-low')
    plt.xlabel("Time (min)")
    plt.ylabel("Reacting Fluid Temperature (°C)")
    plt.savefig('example_6_6_6_response_at_S_low.png')
    plt.savefig('example_6_6_6_response_at_S_low.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch6_ex6/example_6_6_6_response_at_S_low.png')
    plt.legend()
    plt.show()
    return

# calculate the deliverables
if __name__=="__main__":
    deliverables()
