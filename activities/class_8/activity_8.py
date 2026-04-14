"""Generate Concentration Profile for Learning Activity 8"""

# import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc("savefig", dpi=300)

# kinetics parameters from First CoRE Activity 30
dH = 14e3 # cal/mol
k0 = 4.2e15 # cm^3/mol/min
E = 18e3 # cal/mol
T0 = 300 # K
CA0 = 2e-3 # mol/cm^3
CZ0 = 0 # mol/l
V = 1000 # cm^3
Cp = 1.3 # cal/cm^3/K
tf = 80 # min
R = 1.987 # cal/mol/K

# bstr model
def iso_bstr_model_variables():
    # set the initial values
    ind_0 = 0.0
    dep_0 = [CA0*V, CZ0*V]

    # define the stopping criterion
    f_var = 0
    f_val = tf

    # solve the IVODE design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                ,iso_bstr_model_derivatives, odes_are_stiff=True)
    
    # check that a solution was found
    if not(success):
        print(f"!! BSTR error: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]

    # return the profiles
    return t, nA, nZ

def iso_bstr_model_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nZ = dep[1]

    # calculate the reaction rate
    k = k0*np.exp(-E/(R*T0))
    r = k*nA/V*nZ/V + 0.01*k*nA**2/V**2

    # calculate the derivatives
    dnA_dt = -r*V
    dnZ_dt = r*V

    # return the derivatives as a list
    return [dnA_dt, dnZ_dt]

# adiabatic bstr model
def adiab_bstr_variable():
    # set the initial values
    ind_0 = 0.0
    dep_0 = [CA0*V, CZ0*V, T0]

    # define the stopping criterion
    f_var = 0
    f_val = tf

    # solve the IVODE design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                ,adiab_bstr_model_derivatives, odes_are_stiff=True)
    
    # check that a solution was found
    if not(success):
        print(f"!! BSTR error: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]

    # return the profiles
    return t, nA, nZ, T

# adiabatic bstr model derivatives
def adiab_bstr_model_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]

    # calculate the reaction rate
    k = k0*np.exp(-E/(R*T))
    r = k*nA/V*nZ/V + 0.01*k*nA**2/V**2

    # calculate the derivatives
    dnA_dt = -r*V
    dnZ_dt = r*V
    dT_dt = -dH*r/(Cp)

    # return the derivatives as a list
    return [dnA_dt, dnZ_dt, dT_dt]

# deliverables function
def deliverables():
    # get the profiles for the isothermal bstr model
    t_iso, nA_iso, nZ_iso = iso_bstr_model_variables()

    # get the profiles for the adiabatic bstr model
    t_adiab, nA_adiab, nZ_adiab, T_adiab = adiab_bstr_variable()

    # plot the concentration profiles
    plt.figure(1)
    plt.plot(t_iso, nA_iso/V, label="Isothermal")
    plt.plot(t_adiab, nA_adiab/V, label="Adiabatic")
    plt.xlabel("Time")
    plt.ylabel("Concentration of A")
    plt.yticks([])
    plt.xticks([])
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.ylim
    plt.legend()
    plt.savefig("activity_8_concentration_profiles.pdf")
    plt.show()

# execution command
if __name__ == "__main__":
    deliverables()
