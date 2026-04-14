"""Graph for Example 5.5.1 of REB, The Book"""

# import libraries
import numpy as np
from reb_utils import solve_ivodes
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi=300)

# constants from example 9.6.4 of the previous edition of Reaction Engineering Basics
k0 = 2.59E9 # /min
E = 17500. # cal /mol
dH = 8000. # cal /mol
CA_0 = 2. # mol /L
T_0 = 93  + 273.15 # K
Cp = 440. # cal /L /K
V = 4.0 # L

# constant added for this example
tf = 500 # min

# known
Re = 1.987 # cal /mol /K

# calculated
nA_0 = CA_0*V

# adiabatic BSTR model function
def adiab_bstr_model_variables():
    # set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, 0.0, T_0])

	# define the stopping criterion
    f_var = 0
    f_val = tf
    
	# solve the IVODEs
    odes_are_stiff = False
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                    , adiab_bstr_derivatives, odes_are_stiff)
    
    # check that a solution was found
    if not(success):
        print(f"!! BSTR error: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]

    # return all profiles
    return t, nA, nZ, T

# adiabatic BSTR derivatives function
def adiab_bstr_derivatives(ind, dep):
    # extract the individual dependent variables
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]

    # calculate the rate
    CA = nA/V
    k = k0*np.exp(-E/Re/T)
    r = k*CA

    dnAdt = -V*r
    dnZdt = V*r
    dTdt = (-V*r*dH)/V/Cp

    # return the derivatives
    return [dnAdt, dnZdt, dTdt]

# isothermal BSTR model function
def iso_bstr_model_variables():
    # set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, 0.0])

	# define the stopping criterion
    f_var = 0
    f_val = tf
    
	# solve the IVODEs
    odes_are_stiff = False
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                    , iso_bstr_derivatives, odes_are_stiff)
    
    # check that a solution was found
    if not(success):
        print(f"!! BSTR error: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]

    # return all profiles
    return t, nA, nZ

# isothermal BSTR derivatives function
def iso_bstr_derivatives(ind, dep):
    # extract the individual dependent variables
    nA = dep[0]
    nZ = dep[1]

    # calculate the rate
    CA = nA/V
    k = k0*np.exp(-E/Re/T_0)
    r = k*CA

    dnAdt = -V*r
    dnZdt = V*r

    # return the derivatives
    return [dnAdt, dnZdt]

# deliverables function
def deliverables():
    # solve the BSTR design equations
    t_adiab, nA_adiab, nZ_adiab, T_adiab = adiab_bstr_model_variables()
    t_iso, nA_iso, nZ_iso = iso_bstr_model_variables()
    CA_adiab = nA_adiab/V
    CA_iso = nA_iso/V

    # graph the concentration profiles
    plt.figure(1)
    plt.plot(t_adiab, CA_adiab, label = "Adiabatic")
    plt.plot(t_iso, CA_iso, label = "Isothermal")
    plt.xlabel("Time")
    plt.ylabel("Concentration of A")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xticks([])
    plt.yticks([])
    plt.legend()
    plt.savefig("practice_8_full_conc_profile.pdf")
    plt.show(block = False)

    # graph the temperature profiles
    plt.figure(2)
    plt.plot(t_adiab, T_adiab, label = "Adiabatic")
    plt.axhline(T_0, label = "Isothermal", color = "orange")
    plt.xlabel("Time")
    plt.ylabel("Temperature")
    plt.xlim(left=0)
    plt.xticks([])
    plt.yticks([])
    plt.legend()
    plt.savefig("practice_8_full_temp_profile.pdf")
    plt.show()

# execution command
if __name__=="__main__":
    deliverables()
