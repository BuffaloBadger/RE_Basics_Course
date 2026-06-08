"""Calculations for the Class 22 Practice Assignment from REB, The Course"""

# import libraries
import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes
from reb_utils import solve_ates

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
P_in = 45 # psi
Vdot_in = 120 # ft^3 /h
CA_in = 0.025 # lbmol / ft^3
T_in = 120 + 459.67 # R
Tex_in = 75 + 459.67 # R
mDot = 2000 # lb / h
D = 1/12 # ft
L = 125 # ft
U = 150 # BTU /ft^2 / h / R
fD = 0.018
Cp = 8 # BTU /lbmol / R
mu = 1 # lb / ft / h
rho = 57 # lb / ft^3
Cp_ex = 1 # BTU / lb / R
dH_1 = -30500 # BTU / lbmol
k_1_120 = 0.059*3600 # 1 / h
E_1 = 14000 # BTU / lbmol
# known
R = 1.987 # BTU / lbmol / R
# calculated
nDotA_in = CA_in * Vdot_in
Vdot = Vdot_in
k0_1 = k_1_120 * np.exp(E_1 / (R * (120 + 459.67)))
G = 4*Vdot_in*rho / (np.pi * D**2)

# global variables for the current value of the exchange fluid temperature
g_Tex = float('nan')

# PFR reactor function
def pfr_model_variables(Tex):
    # make Tex available to the derivatives function
    global g_Tex
    g_Tex = Tex

    # define the initial values
    ind_0 = 0
    dep_0 = np.array([nDotA_in, 0, T_in, P_in])

    # define the stopping criterion
    f_var = 0
    f_val = L

    # solve the design equations
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , pfr_derivatives, odes_are_stiff=True)
    
    # check for solver issues
    if not success:
        print('')
        print(f'PFR model issue: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit.')

    # return the pfr model variables
    return z, dep[0,:], dep[1,:], dep[2,:], dep[3,:]

# PFR derivatives function
def pfr_derivatives(z, dep):
    # extract the dependent variables
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]
    P = dep[3]

    # calculate the additional unknowns
    k_1 = k0_1 * np.exp(-E_1 / (R * T))
    CA = nA / Vdot
    r_1 = k_1 * CA

    # evaluate the derivatives
    dnA_dz = np.pi*D**2/4*(-r_1)
    dnB_dz = np.pi*D**2/4*(r_1)
    dT_dz = (np.pi*D*U*(g_Tex - T) - np.pi*D**2/4*r_1*dH_1) / ((nA + nZ)*Cp)
    dP_dz = -fD*G**2/(2*D*rho)/32.174/(3600**2)/144 # psi/ft

    # return the derivatives
    return [dnA_dz, dnB_dz, dT_dz, dP_dz]

# coupled unknown residual function
def coupled_unknown_residual(guess):
    # solve the PFR design equations
    z, nDotA, nDotZ, T, P = pfr_model_variables(guess)

    # evaluate and return the residual
    Qdot = np.pi*D*U*(sp.trapezoid(guess - T, z))
    epsilon = Qdot + mDot*Cp_ex*(guess - Tex_in)

    return epsilon

# deliverables function
def deliverables():
    # guess the exchange fluid temperature
    Tex_guess = Tex_in + 5

    # calculate the exchange fluid temperature
    soln, success, message = solve_ates(coupled_unknown_residual, Tex_guess)
    Tex = soln[0]

    # check for solver issues
    if not success:
        print('')
        print(f'Coupled Unknown issue: {message}')
        print('')
        input('Press return to continue of CTRL-C to exit')

    print(f'The outlet exchange fluid temperature is {Tex - 459.6:.0f} °F')
    
    # solve the PFR design equations
    z, nDotA, nDotZ, T, P = pfr_model_variables(Tex)

    # calculate the conversion
    fA = 100*(nDotA_in - nDotA)/nDotA_in

    # generate the requested graphs
    plt.figure(1)
    plt.plot(z, fA)
    plt.xlabel('Axial position, z (ft)')
    plt.xlim(left=0)
    plt.ylabel('Conversion (%)')
    plt.ylim(bottom=0)
    plt.savefig('practice_22_fA_vs_z.pdf')
    plt.show(block=False)

    plt.figure(2)
    plt.plot(z, T-459.6)
    plt.xlabel('Axial position, z (ft)')
    plt.xlim(left=0)
    plt.ylabel('Temperature (°F)')
    plt.savefig('practice_22_T_vs_z.pdf')
    plt.show(block=False)
    
    plt.figure(3)
    plt.plot(z, P)
    plt.xlabel('Axial position, z (ft)')
    plt.xlim(left=0)
    plt.ylabel('Pressure (psi)')
    plt.savefig('practice_22_P_vs_z.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    