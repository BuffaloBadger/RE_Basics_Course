"""Generate graph for Example 5.5.2 of REB, The Book"""

# import libraries
import numpy as np
from reb_utils import solve_ivodes
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi=300)

# global constants
k1_0 = 5.35e10 # 1/min
E1overR = 9000 # K
k2_0 = 4.61e17 # 1/min
E2overR = 15000 # K
T = 350 # K
CA0 = 1 # mol/L
tf = 30 # min

# BSTR model function
def bstr_model_variables():
    # set the initial values
    ind_0 = 0.0
    dep_0 = np.array([CA0, 0.0, 0.0])

	# define the stopping criterion
    f_var = 0
    f_val = tf
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , bstr_derivatives, False)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # return the reactor_model_variables
    return t, dep[0,:], dep[1,:], dep[2,:]


# BSTR derivatives function
def bstr_derivatives(ind, dep):
    # extract the individual dependent variables
    CA = dep[0]
    CD = dep[1]
    CU = dep[2]

    # calculate the rates
    r1 = k1_0 * np.exp(-E1overR / T) * CA
    r2 = k2_0 * np.exp(-E2overR / T) * CD

    # calculate the derivatives
    dCA_dt = -r1
    dCD_dt = r1 - r2
    dCU_dt = r2

    # return the derivatives
    return dCA_dt, dCD_dt, dCU_dt

# deliverables function
def deliverables():
    # solve the BSTR design equations
    t, CA, CD, CU = bstr_model_variables()
    r1 = k1_0 * np.exp(-E1overR / T) * CA
    r2 = k2_0 * np.exp(-E2overR / T) * CD

    # display and save the concentration profiles
    plt.figure(1)
    plt.plot(t, CA, label="A")
    plt.plot(t, CD, label="D")
    plt.plot(t, CU, label="U")
    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.xticks([])
    plt.yticks([])
    plt.ylim(bottom=0)
    plt.legend()
    plt.savefig("example_5_5_2_concentration_profiles.png")
    plt.savefig("../../../RE_Basics/solutions/ch5_ex2/example_5_5_2_concentration_profiles.png")
    plt.show()

    plt.figure(2)
    plt.plot(t, r1, label="Reaction 1")
    plt.plot(t, r2, label="Reaction 2")
    plt.xlabel("Time")
    plt.ylabel("Rate")
    plt.xticks([])
    plt.yticks([])
    plt.ylim(bottom=0)
    plt.legend()
    plt.savefig("example_5_5_2_rate_profiles.png")
    plt.savefig("../../../RE_Basics/solutions/ch5_ex2/example_5_5_2_rate_profiles.png")
    plt.show()

# execution command
if __name__ == "__main__":
    deliverables()
