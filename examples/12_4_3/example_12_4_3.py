"""Calculations for Example 12.4.3 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ates

# constants available to all functions
# given
T0 = 30 + 273.15 # K
CA_0 = 1.0 # mol /l
CB_0 = 1.2 # mol /l
nDotY_0 = 0.0
nDotZ_0 = 0.0
Vdot = 75 # l /min
k0 = 8.72E5 # l /mol /min
E = 7200 # cal /mol
dH = -10700 # cal /mol
Cp = 1.0 # cal /g /K
rho = 1.0E3 # g /l
fA2 = 0.9
# known
R = 1.987 # cal /mol /K
# calculated
nA0 = CA_0*Vdot
nB0 = CB_0*Vdot

# global variables for quantities that can't be passed to the residuals function
global g_nAin, g_nBin, g_nYin, g_nZin, g_Tin, g_fA
g_nAin = float('nan')
g_nBin = float('nan')
g_nYin = float('nan')
g_nZin = float('nan')
g_Tin = float('nan')
g_fA = float('nan')

# deliverables function
def deliverables():
    # choose a range of values for fA1
    fA1_range = np.linspace(0.05,0.85,100)

    # allocate storage for the corresponding volumes
    V1 = np.ones_like(fA1_range)*float('nan')
    Vtot = np.ones_like(fA1_range)*float('nan')

    # define an initial guess for the first fA1 value
    V1guess = 0.01
    V2guess = 100

    # loop through the R1 conversion values
    for i, fA1 in enumerate(fA1_range):
        # solve the design equations for R1
        V1[i], nB1, nY1, nZ1, T1 = cstr_model_variables(nA0, nB0, 0.0, 0.0, T0
            , fA1, V1guess)

        # solve the design equations for R2
        nA1 = nA0*(1-fA1)
        V2, nB2, nY2, nZ2, T2 = cstr_model_variables(nA1, nB1, nY1, nZ1, T1
            , fA2, V2guess)

        # save the total volume
        Vtot[i] = V1[i] + V2

        # use the results as the next guess
        V1guess = V1[i]
        V2guess = V2

    # plot, show, and save the total volume as a function of the R1 conversion
    plt.figure()
    plt.plot(100*fA1_range, Vtot, color='k', ls='-')
    plt.xlabel('Conversion in R1 (%)')
    plt.ylabel('Total Volume (L)')
    plt.savefig('V_vs_fR1.png', dpi=300)
    plt.show()

    # find the minimum total volume and corresponding conversion
    Vmin = np.min(Vtot)
    iMin = np.argmin(Vtot)
    fR1atMin = fA1_range[iMin]
    
    # find the corresponding R1 and R2 volumes
    V1atMin = V1[iMin]
    V2atMin = Vmin - V1atMin

    # tabulate, show and save the results
    item = ['Minimum Total Volume', 'R1 Volume', 'R2 Volume', 'R1 Conversion']
    value = [Vmin, V1atMin, V2atMin, fR1atMin]
    units = ['L', 'L', 'L', '%']
    resultsDf = pd.DataFrame({'Item':item, 'Value':value, 'Units':units})
    print('')
    print(resultsDf)
    print('')
    resultsDf.to_csv('results.csv', index=False)

# CSTR model function
def cstr_model_variables(nAin, nBin, nYin, nZin, Tin, fA, Vguess):
    # make the inlet molar flow rates, inlet temperature, and conversion 
    # available to the residuals function
    global g_nAin, g_nBin, g_nYin, g_nZin, g_Tin, g_fA
    g_nAin = nAin
    g_nBin = nBin
    g_nYin = nYin
    g_nZin = nZin
    g_Tin = Tin
    g_fA = fA

    # define initial guesses
    nBguess = nBin
    nYguess = nYin
    nZguess = nZin
    Tguess = Tin + 5
    guess = np.array([Vguess, nBguess, nYguess, nZguess, Tguess])
     
	# solve the ATEs
    soln, success, message = solve_ates(cstr_residuals, guess)

    # check for solver issues
    if not(success):
        print('')
        print(f'Second CSTR solver issue: {message}')
        print('')
        input('Press return to continue of CTRL-C to exit')

    # extract the results
    V = soln[0]
    nB = soln[1]
    nY = soln[2]
    nZ = soln[3]
    T = soln[4]

    # return the solution
    return V, nB, nY, nZ, T

# CSTR residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    V = guess[0]
    nBout = guess[1]
    nYout = guess[2]
    nZout = guess[3]
    Tout = guess[4]

    # calculate additional unknowns
    k = k0*np.exp(-E/R/Tout)
    nAout = nA0*(1-g_fA)
    CA = nAout/Vdot
    CB = nBout/Vdot
    r = k*CA*CB

    # evaluate the residuals
    eps1 = g_nAin - nAout - V*r
    eps2 = g_nBin - nBout - V*r
    eps3 = g_nYin - nYout + V*r
    eps4 = g_nZin - nZout + V*r
    eps5 = rho*Vdot*Cp*(Tout - g_Tin) + V*r*dH

    # return the residuals as an array
    return np.array([eps1, eps2, eps3, eps4, eps5])

if __name__=="__main__":
    deliverables()
