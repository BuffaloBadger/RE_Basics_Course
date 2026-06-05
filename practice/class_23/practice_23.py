"""Calculations for the Class 23 Practice Assignment from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
DeA = 2.7e-7 # cm^2 /s
k_1 = 0.019# /s
Vdot_in = 1000/60 # cm^3 /s
CA_in = 1E-3 # mol /cm^3
fA = 0.9
Dp_values = np.array([0.04, 0.27, 0.55, 0.77])/10 # cm
# known
# calculated
nDotA_in = Vdot_in*CA_in
nDotA_out = nDotA_in*(1-fA)
Vdot = Vdot_in

# global variable for the current particle diameter
g_Dp = float('nan')

# PFR reactor function
def pfr_model_variables(Dp):
    # make the particle diameter available to the derivatives function
    global g_Dp
    g_Dp = Dp

    # define the initial values
    ind_0 = 0
    dep_0 = np.array([nDotA_in])

    # define the stopping criterion
    f_var = 1
    f_val = nDotA_out

    # solve the PFR design equations
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            ,pfr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f'PFR model function issue: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit.')
    
    # return the pfr model variables
    return z, dep[0,:]

# PFR derivatives function
def pfr_derivatives(ind, dep):
    # extract the dependendent variable
    nDotA = dep[0]

    # calculate the additional unknowns
    CA = nDotA/Vdot
    phi = g_Dp/2*np.sqrt(k_1/DeA)
    effectiveness = 3/phi * (1/np.tanh(phi) - 1/phi)
    r = k_1*CA

    # evaluate and return the derivative
    return -effectiveness*r

# deliverables function
def deliverables():
    # allocate storage for the PFR volumes
    Vpfr = np.ones_like(Dp_values)*float('nan')

    # loop through the Dp values
    for i, Dp in enumerate(Dp_values):
        # solve the design equations
        V, nDotA = pfr_model_variables(Dp)

        # save the PFR volume
        Vpfr[i] = V[-1]

    # tabulate, show, and save the results
    results_df = pd.DataFrame({'Particle Diameter (mm)':10*Dp_values
                               ,'PFR Bed Volume (cc)':Vpfr})
    print('')
    print(results_df)
    print('')
    results_df.to_csv('practice_23_results.csv',index=False)

# execution command
if __name__ == '__main__':
    deliverables()
    