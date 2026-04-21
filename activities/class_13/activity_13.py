"""Calculations for the Class 13 Learning Activity in REB, The Course"""

# import libraries
import numpy as np
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all function
# given
V = 4.430e3 # cm3
CS_in = 0.04 # g/cm3
CX_in = 0 # g/cm3
Vmax = 0.014 # /min
Km = 0.001 # g/cm3
Vdot_in = 61 # cm3/min
CS_0 = 2.97e-2 # g/cm3
CX_0 = 4.68e-3 # g/cm3
# calculated
Vdot = Vdot_in

# cstr model function
def cstr_model_variables():
    # set the initial values
    ind_0 = 0
    dep_0 = np.array([CS_0, CX_0]) *Vdot_in

    # set the stopping criterion
    f_var = 0
    f_val = 40000

    # solve the cstr design equations
    odes_are_stiff = False
    ind, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                , cstr_derivatives, odes_are_stiff)
    
    # check for solver issues
    if not success:
        print('')
        print(f"CSTR model function issue: {message}")
        print('')
        input('Press return to continue.')
    
    # return the cstr model variables
    return ind, dep[0,:], dep[1,:] # t, mS, mX

# cstr derivatives function
def cstr_derivatives(ind, dep):
    # extract the dependent variables
    mS = dep[0]
    mX = dep[1]

    # calculate the additional unknowns
    CS = mS/Vdot
    CX = mX/Vdot
    r1 = Vmax*CS*CX/(Km + CS)
    mS_in = CS_in*Vdot_in
    mX_in = CX_in*Vdot_in

    # evaluate the derivatives
    dmSdt = Vdot_in/V*(mS_in - mS -2.2*V*r1)
    dmXdt = Vdot/V*(mX_in - mX + V*r1)

    # return the derivatives
    return dmSdt, dmXdt

# deliverables function
def deliverables():
    # solve the cstr design equations
    t, mS, mX = cstr_model_variables()

    # generate, show and save a plot of CX vs time
    CX = mX/Vdot
    plt.figure(1)
    plt.plot(t, CX)
    plt.xlabel("Time (min)")
    plt.ylabel("Outlet Cell Mass Concentration (g cm$^{-3}$)")
    plt.ylim(bottom=0)
    plt.legend()
    plt.savefig("activity_13_cell_mass_vs_t.png")
    plt.savefig("activity_13_cell_mass_vs_t.pdf")
    plt.show()

# execution command
if __name__=="__main__":
    deliverables()
