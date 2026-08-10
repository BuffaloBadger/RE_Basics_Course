"""Calculations for Example 12.4.1 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ivodes
from reb_utils import solve_ates

# constants available to all functions
# given
Vcstr = 350.0 # L
Vpfr = 350.0 # L
CAfeed = 2.5 # mol /L
VdotFeed = 100.0 # L /min
Tfeed = 38 + 273.15 # K
dH1 = -21500 # cal /mol
dH2 = -24000 # cal /mol
Cp = 1.0E3 # cal /L /K
k0_1 = 1.2E5 # /min
E_1 = 9100 # cal/mol
k0_2 = 2.17E7 # L /mol /min
E_2 = 13400 # cal /mol
# known
R = 1.987 # cal /mol
# calculated
nDotAfeed = VdotFeed*CAfeed

# global variables to make the current inlet molar flow rates and temperature 
# available to the CSTR residuals function
g_nAin = float('nan')
g_nDin = float('nan')
g_nUin = float('nan')
gTin = float('nan')

# CSTR model function
def cstr_model_variables(nA, nD, nU, T):
    # make the inlet quantities available to the residuals function
    global g_nAin, g_nDin, g_nUin, gTin
    g_nAin = nA
    g_nDin = nD
    g_nUin = nU
    gTin = T

    # define an initial guess
    initGuess = np.array([nA, nD, nU, T + 10])

    # solve the CSTR design equations
    soln, success, message = solve_ates(cstr_residuals, initGuess)

    # check for solver issues
    if not success:
        print('')
        print(f'CSTR model function solver issue: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit.')

    # return the CSTR model variables
    return soln[0], soln[1], soln[2], soln[3]

# CSTR residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nDotAout = guess[0]
    nDotDout = guess[1]
    nDotUout = guess[2]
    Tout = guess[3]

    # calculate the additional unknowns
    k1 = k0_1*np.exp(-E_1/R/Tout)
    k2 = k0_2*np.exp(-E_2/R/Tout)
    CA = nDotAout/VdotFeed
    r1 = k1*CA
    r2 = k2*CA**2

    # evaluate the residuals
    epsilon_1 = g_nAin - nDotAout + (-r1 -r2)*Vcstr
    epsilon_2 = g_nDin - nDotDout + r1*Vcstr
    epsilon_3 = g_nUin - nDotUout + r2*Vcstr
    epsilon_4 = -VdotFeed*Cp*(Tout - gTin) - Vcstr*(r1*dH1 + r2*dH2)

    # return the residuals as an array
    return np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4])

# PFR model function
def pfr_model_variables(nA, nD, nU, T):
    # define the initial values
    ind0 = 0.0
    dep0 = np.array([nA, nD, nU, T])

    # define the stopping criterion
    fVar = 0
    fVal = Vpfr

    # solve the PFR design equations
    V, dep, success, message = solve_ivodes(ind0, dep0, fVar, fVal
        , pfr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not success:
        print('')
        print(f'PFR model function solver issue: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit.')

    # return the PFR model variables
    return V, dep[0,:], dep[1,:], dep[2,:], dep[3,:]

# PFRderivatives function
def pfr_derivatives(ind, dep):
	# extract the dependent variables
    nDotA = dep[0]
    nDotD = dep[1]
    nDotU = dep[2]
    T = dep[3]

	# calculate the additional unknowns
    k1 = k0_1*np.exp(-E_1/R/T)
    k2 = k0_2*np.exp(-E_2/R/T)
    CA = nDotA/VdotFeed
    r1 = k1*CA
    r2 = k2*CA**2

	# evaluate the derivatives
    dnAdV = -r1 -r2
    dnDdV = r1
    dnUdV = r2
    dTdV = -(r1*dH1 + r2*dH2)/VdotFeed/Cp

	# return the derivatives
    return dnAdV, dnDdV, dnUdV, dTdV

# deliverables
def deliverables():
    # case a - solve the CSTR design equations
    nAoutCSTR, nDoutCSTR, nUoutCSTR, ToutCSTR = cstr_model_variables(nDotAfeed
        ,0.0, 0.0, Tfeed)

    # solve the PFR design equations
    V, nDotA, nDotD, nDotU, T = pfr_model_variables(nAoutCSTR, nDoutCSTR
        , nUoutCSTR, ToutCSTR)

    # calculate the quantities of interest
    fAcase_a = 100.0*(nDotAfeed - nDotA[-1])/nDotAfeed
    selCase_a = nDotD[-1]/nDotU[-1]
    Ta = T[-1] - 273.15

    # case b - solve the PFR design equations
    V, nDotA, nDotD, nDotU, T = pfr_model_variables(nDotAfeed, 0.0, 0.0, Tfeed)

    # solve the CSTR design equations
    nAoutCSTR, nDoutCSTR, nUoutCSTR, ToutCSTR = cstr_model_variables(nDotA[-1]
        , nDotD[-1], nDotU[-1], T[-1])

    # calculate the other quantities of interest
    fAcase_b = 100.0*(nDotAfeed - nAoutCSTR)/nDotAfeed
    selCase_b = nDoutCSTR/nUoutCSTR
    Tb = ToutCSTR - 273.15

    # tabulate, show, and save the results
    results = [['CSTR then PFR', fAcase_a, selCase_a, Ta]
               ,['PFR then CSTR', fAcase_b, selCase_b, Tb]]
    results_df = pd.DataFrame(results, columns=['Reactors','conversion'
        ,'selectivity', 'temperature'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('results.csv',index=False)

    return

# execution command
if __name__=="__main__":
    deliverables()
