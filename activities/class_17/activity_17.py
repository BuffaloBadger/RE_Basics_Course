"""Calculations for the Class 17 Learning Activity in REB, The Course"""

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
k01 = 265 # L /mol /min
E1 = 73000 # J /mol
dH298 = -165000 # J /mol
P0 = 1000 # Torr
T0 = 1225 # K
# Constants from heat capacity expressions Cpi = ai + bi*T
aA = 28
aY = 26
aZ = 30
bA = 0.05
bY = 0.01
bZ = 0.005
V = 1 # basis
Tf = np.array([1235,1325]) # K
# known
Re = 8.314 # J /mol /K
Rw = 62.36367 # L Torr /mol \K
# calculated
nA0 = P0*V/(Rw*T0)

# BSTR reactor function
def bstr_model_variables(Tf):
    # define the initial values
    ind_0 = 0
    dep_0 = np.array([nA0, 0, 0, T0, P0])

    # define the stopping criterion
    f_var = 4
    f_val = Tf

    # solve the design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            ,bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model function issue: {message}")
        print('')
        input('Press return to continue.')

    # return the bstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:]

# BSTR derivatives function
def bstr_derivatives(t, dep):
    # extract the dependent variables
    nA = dep[0]
    nY = dep[1]
    nZ = dep[2]
    T = dep[3]
    P = dep[4]

    # calculate the additional unknowns
    Cp_A = aA + bA*T
    Cp_Y = aY + bY*T
    Cp_Z = aZ + bZ*T
    dH_1 = dH298 + (aZ + 2*aY - 2*aA)*(T-298) + 0.5*(bZ + 2*bY - 2*bA)*(T**2-298**2)
    r_1 = k01*np.exp(-E1/Re/T)*(nA/V)**2

	# Create mass matrix, setting all elements to zero
    mass_matrix = np.zeros((5,5))

    # Add 1 on the diagonal for the first 3 rows
    mass_matrix[0,0] = 1.0
    mass_matrix[1,1] = 1.0
    mass_matrix[2,2] = 1.0

    # Add the elements for the energy balance
    mass_matrix[3,3] = nA*Cp_A + nY*Cp_Y + nZ*Cp_Z
    mass_matrix[3,4] = -V*Re/Rw

    # Add the elements for the ideal gas law equation
    mass_matrix[4,0] = Rw*T
    mass_matrix[4,1] = Rw*T
    mass_matrix[4,2] = Rw*T
    mass_matrix[4,3] = Rw*(nA + nY + nZ)
    mass_matrix[4,4] = -V

    # Create right side vector
    rhs1 = -2*V*r_1
    rhs2 = 2*V*r_1
    rhs3 = V*r_1
    rhs4 = -V*r_1*dH_1
    rhs5 = 0.0
    rhs = np.array([rhs1, rhs2, rhs3, rhs4, rhs5])

    # Evaluate the derivatives
    derivs = sp.linalg.solve(mass_matrix, rhs)

    # Return the derivatives
    return derivs

# deliverables function
def deliverables():
    # allocate storage for the deliverables
    tf = np.ones_like(Tf) * float('NaN')
    fA = np.ones_like(Tf) * float('NaN')

    # loop through the Tf values
    for n, Tfn in enumerate(Tf):
        # solve the bstr design equations
        t, nA, nY, nZ, T, P = bstr_model_variables(Tfn)

        # save the final time and the conversion
        tf[n] = t[-1]
        fA[n] = 100*(nA0 - nA[-1])/nA0
    
    # tabulate, show, and save the results
    results_df = pd.DataFrame({"Final T (K)" : Tf, "Time (min)" : tf, "Conversion (%)" : fA})
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_17_results.csv',index=False)

# execution command
if __name__ == '__main__':
    deliverables()
    