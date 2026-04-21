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
Cp = 0.35*4.184 # J /g /K
rho = 0.93 # g /cm^3
dH = -20000.0 # J /mol
k0 = 3.24E12 # cm^3 /mol /s
E = 105000.0 # J /mol
# known
R = 8.314 # J /mol /K
# calculated
Vdot = Vdot_in
nA_in = CA_in*Vdot_in
nB_in = CB_in*Vdot_in

# cstr model function
def cstr_model_variables(T_0, t_f):
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([0, 0, 0, 0, T_0])

	# define the stopping criterion
    f_var = 0
    f_val = t_f
     
	# solve the IVODEs
    odes_are_stiff = False
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , cstr_derivatives, odes_are_stiff)

    # check that a solution was found
    if not(success):
        print('')
        print(f"CSTR model function issue: {message}")
        print('')
        input('Press enter to continue.')

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

	# calculate the additional unknowns
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

    # solve the cstr design equations for an initial temperature of 180 °V
    [t_180, nA_180, nB, nY, nZ, T_180] = cstr_model_variables(180 + 273.15, t_f)

    # solve the cstr design equations for an initial temperature of 190 °V
    [t_190, nA_190, nB, nY, nZ, T_190] = cstr_model_variables(190 + 273.15, t_f)

    # create, display and save the graphs
    plt.figure(1) 
    plt.plot(t_180,T_180 - 273.15, label = '180 °C')
    plt.plot(t_190,T_190 - 273.15, label = '190 °C')
    plt.xlabel("Time (min)")
    plt.ylabel("Reacting Fluid Temperature (°C)")
    plt.legend()
    plt.savefig('activity_12_T_profiles.png')
    plt.savefig('activity_12_T_profiles.pdf')
    plt.show(block = False)

    CA_180 = nA_180/Vdot
    CA_190 = nA_190/Vdot
    plt.figure(2)
    plt.plot(t_180, CA_180, label = '180 °C')
    plt.plot(t_190, CA_190, label = '190 °C')
    plt.xlabel("Time (min)")
    plt.ylabel("Concentration of A (mol cm$^{-3}$)")
    plt.legend()
    plt.savefig('activity_12_CA_profiles.png')
    plt.savefig('activity_12_CA_profiles.pdf')
    plt.show()
    return

# calculate the deliverables
if __name__=="__main__":
    deliverables()
