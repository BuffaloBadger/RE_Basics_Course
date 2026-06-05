"""Calculations for the Class 22 Practice Assignment from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
y_D_in = 0.08
y_O_in = 0.13
y_I_in = 0.79
P_in = 1.0 # atm
T_in = 370.0 # °C
T_in = (T_in + 273.15) * 1.8 # °R
nDot_in = 0.149 # lbmol /s
nDot_in = nDot_in * 3600 # lbmol /h
mu = 0.09 # lb /h /ft
Dp = 0.25 # in
Dp = Dp / 12 # ft
eps = 0.4
rho_bed = 0.6 # g /cm^3
rho_bed = rho_bed / 453.6 *28317 # lb /ft^3
D = 6.0 # ft
f_D = 0.81 
k0_f = 1.745E5 # mol /s /g_cat /atm^1.5
k0_f = k0_f * 3600 # lbmol /h lb_cat /atm^1.5
E_f = 31000 # cal /mol
E_f = E_f   * 0.00397*453.6 # BTU /lbmol
k0_r = 7.59E9 # mol /s /g_cat /atm
k0_r = k0_r * 3600 # lbmol /h lb_cat /atm
E_r = 53600 # cal /mol
E_r = E_r * 0.00397*453.6 # BTU /lbmol
hf298 = np.array([-70950, 0, -94470, 0]) # cal /mol
hf298 = hf298 * 0.00397*453.6 # BTU /lbmol
alpha = np.array([5.697, 6.713, 12.13, 7.44]) # cal /mol /K
alpha = alpha * 0.00397*453.6/1.8 # BTU /lbmol /°R
beta = np.array([0.016, -8.790E-07, 0.00812, -0.00324]) # cal /mol /K**2
beta = beta * 0.00397*453.6/1.8**2 # BTU /lbmol /°R**2
gamma = np.array([-1.185E-05, 4.175E-06, 0, 6.4E-06]) # cal /mol /K**3
gamma = gamma * 0.00397*453.6/1.8**3 # BTU /lbmol /°R**3
delta = np.array([3.172E-09, -2.544E-09, 0.0, -2.790E-09]) # cal /mol /K**4
delta = delta * 0.00397*453.6/1.8**4 # BTU /lbmol /°R**4
# known
mw_D = 64 # lbm/lbmol
mw_O = 32 # lbm/lbmol
mw_I = 28 # lbm/lbmol
Ren = 1.986 # BTU /lbmol /°R
Rpv = 0.7302 # ft3 atm / °R / lbmol
# calculated
nDot_D_in = y_D_in*nDot_in
nDot_O_in = y_O_in*nDot_in
nDot_I_in = y_I_in*nDot_in
nDot_D_out = nDot_D_in * (1 - f_D)
G = 4*(nDot_D_in*mw_D + nDot_O_in*mw_O + nDot_I_in*mw_I)/(np.pi*D**2) # lbm /h /ft2
dH_298 = hf298[2] - 0.5*hf298[1] - hf298[0] # BTU /lbmol
dalpha = alpha[2] - 0.5*alpha[1] - alpha[0] # BTU /lbmol /°R
dbeta = beta[2] - 0.5*beta[1] - beta[0] # BTU /lbmol /°R**2
dgamma = gamma[2] - 0.5*gamma[1] - gamma[0] # BTU /lbmol /°R**3
ddelta = delta[2] - 0.5*delta[1] - delta[0] # BTU /lbmol /°R**4

# parameter and global variable for its current value
D_values = np.array([6, 8, 10]) # ft
g_D = float('nan')

# PFR reactor function
def pfr_model_variables(D):
    # make the current value of D available to the derivatives function
    global g_D
    g_D = D

    # define the initial values
    ind_0 = 0
    dep_0 = np.array([nDot_D_in, nDot_O_in, 0, nDot_I_in, T_in, P_in])

    # define the stopping criterion
    f_var = 1
    f_val = nDot_D_out

    # solve the design equations
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            ,pfr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f'PFR model issue: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit.')
    
    # return the pfr model variables
    return z, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:]

# PFR derivatives function
def pfr_derivatives(ind, dep):
    # use the current value of D
    D = g_D

    # extract the dependent variables
    nDot_D = dep[0]
    nDot_O = dep[1]
    nDot_T = dep[2]
    nDot_I = dep[3]
    T = dep[4]
    P = dep[5]

    # calculate the additional unknowns
    nDot = nDot_D + nDot_O + nDot_T + nDot_I
    vDot = nDot*Rpv*T/P
    rho = G * np.pi*D**2/4 / vDot
    dT = T - 536.4
    dT2 = T**2 - 536.4**2
    dT3 = T**3 - 536.4**3
    dT4 = T**4 - 536.4**4
    dH = dH_298 + dalpha*dT + dbeta/2.*dT2 + dgamma/3.*dT3 + ddelta/4.*dT4
    CpD = alpha[0] + beta[0]*T + gamma[0]*T**2 + delta[0]*T**3
    CpO = alpha[1] + beta[1]*T + gamma[1]*T**2 + delta[1]*T**3
    CpT = alpha[2] + beta[2]*T + gamma[2]*T**2 + delta[2]*T**3
    CpI = alpha[3] + beta[3]*T + gamma[3]*T**2 + delta[3]*T**3
    PD = nDot_D*P/nDot
    PO = nDot_O*P/nDot
    PT = nDot_T*P/nDot
    k_f = k0_f*np.exp(-E_f/Ren/T)
    k_r = k0_r*np.exp(-E_r/Ren/T)
    r = rho_bed*(k_f*PD*PO - k_r*PT*np.sqrt(PO))/np.sqrt(PD)

    # evaluate the derivatives
    dnDotDdz = -np.pi*D**2/4*r
    dnDotOdz = -0.5*np.pi*D**2/4*r
    dnDotTdz = np.pi*D**2/4*r
    dnDotIdz = 0.0
    dTdz = (-np.pi*D**2/4*r*dH)/(nDot_D*CpD + nDot_O*CpO
            + nDot_T*CpT + nDot_I*CpI)
    dPdz = 0

    # return the pfr derivatives
    return [dnDotDdz, dnDotOdz, dnDotTdz, dnDotIdz, dTdz, dPdz]

# deliverables function
def deliverables():
    # allocate storage for the quantities of interest
    V_cat = np.ones_like(D_values) * float('nan')
    P_out = np.ones_like(D_values) * float('nan')

    # define a figure for the T vs. z graph
    plt.figure(1)

    # loop through the D values
    for i, D in enumerate(D_values):
        # solve the PFR design equations
        z, nDot_D, nDot_O, nDot_T, nDot_I, T, P = pfr_model_variables(D)

        # calculate the quantities of interest
        V_cat[i] = np.pi*D**2/4*z[-1]
        P_out[i] = P[-1]

    # tabulate, show, and save the results
    results_df = pd.DataFrame({'Diameter (ft)': D_values, 'Catalyst Volume (ft^3)' : V_cat
                            ,'Pressure Drop (atm)': P_in - P_out})
    print('')
    print(results_df)
    print('')
    results_df.to_csv('practice_22_discussion_results.csv',index=False)

# execution command
if __name__ == '__main__':
    deliverables()
