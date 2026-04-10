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
dH = -12200. # cal /mol
CA_0 = 2. # mol /L
T_0 = 23  + 273.15 # K
Cp = 440. # cal /L /K
V = 4.0 # L

# constant added for this example
tf = 1000.0 # min

# known
Re = 1.987 # cal /mol /K

# calculated
nA_0 = CA_0*V

# BSTR model function
def bstr_model_variables():
    # set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, 0.0, T_0])

	# define the stopping criterion
    f_var = 0
    f_val = tf
    
	# solve the IVODEs
    odes_are_stiff = False
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                    , bstr_derivatives, odes_are_stiff)
    
    # check that a solution was found
    if not(success):
        print(f"!! BSTR error: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]

    # return all profiles
    return t, nA, nZ, T

# BSTR derivatives function
def bstr_derivatives(ind, dep):
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

# deliverables function
def deliverables():
    # solve the BSTR design equations
    t, nA, nZ, T = bstr_model_variables()
    CA = nA/V
    r = k0*np.exp(-E/Re/T)*CA

    # graph the concentration profile
    plt.figure(1)
    plt.plot(t, CA)
    plt.xlabel("Time")
    plt.ylabel("Concentration of A")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xticks([])
    plt.yticks([])
    plt.savefig("example_5_5_1_full_conc_profile.png")
    plt.savefig("../../../RE_Basics/solutions/ch5_ex1/example_5_5_1_full_conc_profile.png")
    plt.show(block=False)

    # graph the temperature profile
    plt.figure(2)
    plt.plot(t, T)
    plt.xlabel("Time")
    plt.ylabel("Temperature")
    plt.xlim(left=0)
    plt.xticks([])
    plt.yticks([])
    plt.savefig("example_5_5_1_full_temp_profile.png")
    plt.savefig("../../../RE_Basics/solutions/ch5_ex1/example_5_5_1_full_temp_profile.png")
    plt.show(block=False)

    # graph the rate profile
    plt.figure(3)
    plt.plot(t, r)
    plt.xlabel("Time")
    plt.ylabel("Net Rate of Reaction")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xticks([])
    plt.yticks([])
    plt.savefig("example_5_5_1_full_rate_profile.png")
    plt.savefig("../../../RE_Basics/solutions/ch5_ex1/example_5_5_1_full_rate_profile.png")
    plt.show(block=False)

    # graph the concentration profile only where the rate is decreasing
    i_opt = np.argmax(r)
    plt.figure(4)
    plt.plot(t[i_opt:-1], CA[i_opt:-1])
    plt.xlabel("Time")
    plt.ylabel("Concentration of A")
    plt.xlim(left=t[i_opt])
    plt.ylim(bottom=0)
    plt.xticks([])
    plt.yticks([])
    plt.savefig("example_5_5_1_conc_profile_when_r_decreasing.png")
    plt.savefig("../../../RE_Basics/solutions/ch5_ex1/example_5_5_1_conc_profile_when_r_decreasing.png")
    plt.show(block=False)

    # graph the temperature profile only where the rate is decreasing
    plt.figure(5)
    plt.plot(t[i_opt:-1], T[i_opt:-1])
    plt.xlabel("Time")
    plt.ylabel("Temperature")
    plt.xlim(left=t[i_opt])
    plt.xticks([])
    plt.yticks([])
    plt.savefig("example_5_5_1_temp_profile_when_r_decreasing.png")
    plt.savefig("../../../RE_Basics/solutions/ch5_ex1/example_5_5_1_temp_profile_when_r_decreasing.png")
    plt.show(block=False)

    # graph approach to equilibrium
    plt.figure(6)
    plt.rcParams['axes.spines.left'] = False
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.spines.bottom'] = False
    plt.plot(t[i_opt:-1], CA[i_opt:-1], label="decreasing concave upward")
    plt.plot(t[i_opt:-1], -CA[i_opt:-1], label="increasing concave downward")
    plt.xlabel("Time")
    plt.ylabel("Concentration or Temperature")
    plt.xticks([])
    plt.yticks([])
    plt.savefig("concave_up_and_down.png")
    plt.savefig("../../../RE_Basics/figures/concave_up_and_down.png")
    plt.show()

# execution command
if __name__=="__main__":
    deliverables()
