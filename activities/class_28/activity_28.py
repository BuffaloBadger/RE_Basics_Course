"""Calculations for the Class 28 Learning Activity from REB, The Course"""

import pandas as pd
import numpy as np
from reb_utils import solve_ivodes
from reb_utils import fit_to_SR_data
import matplotlib.pyplot as plt
import os

# global constants available to all functions
# given
P = 1.0 # atm
mTotal = 3.0 # g
VFR_in = 0.85 # L/min
T = 400. + 273.15 # K
K = 12.2
# calculated
nTotalIn = VFR_in/22.4 # mol/min

# experimental data as arrays
df = pd.read_csv('activity_28_data.csv')
yAin = df['yA'].to_numpy()
yBin = df['yB'].to_numpy()
yYin = df['yY'].to_numpy()
yZin = df['yZ'].to_numpy()
PAout = df['PA'].to_numpy()

# global variables for quantities that can't be passed to the derivatives 
# function
global gk, gKA, gKB, gKY, gKZ
gk = float('nan')
gKA = float('nan')
gKB = float('nan')
gKY = float('nan')
gKZ = float('nan')

# PFR model function
def pfr_model_variables(yAin, yBin, yYin, yZin, k, KA, KB, KY, KZ):
    # make the rate expression parameters available to the derivatives function
    global gk, gKA, gKB, gKY, gKZ
    gk = k
    gKA = KA
    gKB = KB
    gKY = KY
    gKZ = KZ

    #initial values and stopping criterion
    ind0 = 0
    nAin = yAin*nTotalIn
    nBin = yBin*nTotalIn
    nYin = yYin*nTotalIn
    nZin = yZin*nTotalIn
    dep0 =[nAin, nBin, nYin, nZin]
    f_var = 0
    f_val = mTotal

    # solve the PFR design equations
    m, dep, success, message = solve_ivodes(ind0, dep0, f_var, f_val
        , pfr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not(success):
        print('')
        print(f"PFR model solver issue: {message}")
        print('')
        input('Press return to continue of CTRL-C to exit.')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]

    # return the profiles
    return m, nA, nB, nY, nZ

# PFR derivatives function
def pfr_derivatives (m, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nY = dep[2]
    nZ = dep[3]

    # calculate the additional unknowns
    ntot = nA + nB + nY + nZ
    PA = nA/ntot*P
    PB = nB/ntot*P
    PY = nY/ntot*P
    PZ = nZ/ntot*P
    r = gk*PA*PB/(1 + gKA*PA + gKB*PB + gKY*PY + gKZ*PZ)\
        *(1 - PY*PZ/K/PA/PB)
    
    # evaluate the PFR derivatives
    dnAdm = -r
    dnBdm = -r
    dnYdm = r
    dnZdm =r

    # return the PFR derivatives as an array
    return np.array([dnAdm, dnBdm, dnYdm, dnZdm])

# predicted responses function
def predicted_responses(adj_inputs, paramsk, paramsA, paramsB, paramsY
        , paramsZ):
    # allocate storage for the responses
    PAout_pred = np.ones_like(PAout) * float('nan')

    # extract the rate expression parameters
    k = 10**paramsk
    KA = 10**paramsA
    KB = 10**paramsB
    KY = 10**paramsY
    KZ = 10**paramsZ

    # loop through the data points
    for i, input in enumerate(adj_inputs):
        # extract the adjusted inputs
        yAin = input[0]
        yBin = input[1]
        yYin = input[2]
        yZin = input[3]

        # solve the PFR design equations
        m, nA, nB, nY, nZ = pfr_model_variables(yAin, yBin, yYin, yZin, k
            , KA, KB, KY, KZ)
        
        # calculate the response
        nAout = nA[-1]
        nBout = nB[-1]
        nYout = nY[-1]
        nZout = nZ[-1]
        PAout_pred[i] = nAout/(nAout + nBout + nYout + nZout)*P
    
    # return the responses
    return PAout_pred

# deliverables function
def deliverables():
    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([yAin, yBin, yYin, yZin]))

    # make a guess the parameters
    #par_guess = [0.0, 0.0, 0.0, 0.0, 0.0] # underpredicted the response
    #par_guess = [-2.0, 0.0, 0.0, 0.0, 0.0] # bad confidence intervals
    par_guess = [0.0, 2.0, 2.0, 2.0, 2.0]

    guessing=False
    if guessing:
        # calculate the predicted responses
        PAout_pred = predicted_responses(adj_inputs, par_guess[0], par_guess[1]
            , par_guess[2], par_guess[3], par_guess[4])
    else:
        # estimate the parameters
        params, params_ci, r_squared, PAout_pred = fit_to_SR_data(par_guess
            , adj_inputs , PAout, predicted_responses, use_rel_error=False)

        # extract the results
        k = 10**params[0]
        k_CI = 10**params_ci[0,:]
        KA = 10**params[1]
        KA_CI = 10**params_ci[1,:]
        KB = 10**params[2]
        KB_CI = 10**params_ci[2,:]
        KY = 10**params[3]
        KY_CI = 10**params_ci[3,:]
        KZ = 10**params[4]
        KZ_CI = 10**params_ci[4,:]
    
        # tabulate, show, and save the results to a .csv file
        kUnits = "mol g^-1^ min^-1^ atm^-2^"
        KUnits = "atm^-1^"
        data = [['k', f'{k:.3g}', f'{k_CI[0]:.3g}', f'{k_CI[1]:.3g}', kUnits],
            ['KA', f'{KA:.2g}', f'{KA_CI[0]:.2g}', f'{KA_CI[1]:.2g}', KUnits],
            ['KB', f'{KB:.2g}', f'{KB_CI[0]:.2g}', f'{KB_CI[1]:.2g}', KUnits],
            ['KY', f'{KY:.2g}', f'{KY_CI[0]:.2g}', f'{KY_CI[1]:.2g}', KUnits],
            ['KZ', f'{KZ:.2g}', f'{KZ_CI[0]:.2g}', f'{KZ_CI[1]:.2g}', KUnits],
            ['R_squared', f'{r_squared:.3g}', '', '', '']]
        results_df = pd.DataFrame(data, columns=['Parameter', 'Value'
                , 'lower limit', 'upper limit', 'Units'])
        print(" ")
        print(results_df)
        results_df.to_csv('results.csv', index=False)
        
    # calculate the experiment residuals
    epsilon_expt = PAout - PAout_pred

    # make sure a results folder and a png folder exist
    if not os.path.isdir('png'):
        # create the folder
        os.makedirs('./png')

    # create, show, and save the parity plot
    plt.figure() 
    plt.plot(PAout, PAout_pred, color = 'k', marker='o', markerfacecolor='none'
        , ls='', label = 'Data')
    plt.plot([min(PAout),max(PAout)],[min(PAout),max(PAout)]
             , color = 'k', ls = '-', label = 'Parity Line')
    plt.xlabel("$P_{A,1}$ (atm)")
    plt.ylabel("$P_{A,1,pred}$ (atm)")
    plt.legend
    plt.tight_layout()
    plt.savefig('png/parity.png', dpi=300)
    plt.show(block=False)

    # create, show, and save residuals plots
    plt.figure() 
    plt.plot(yAin, epsilon_expt, color = 'k', marker='o', markerfacecolor='none'
        , ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$y_{A,0}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('png/A_residual.png', dpi=300)
    plt.show(block=False)

    plt.figure() 
    plt.plot(yBin, epsilon_expt, color = 'k', marker='o', markerfacecolor='none'
        , ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$y_{B,0}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('png/B_residual.png', dpi=300)
    plt.show(block=False)

    plt.figure() 
    plt.plot(yYin, epsilon_expt, color = 'k', marker='o', markerfacecolor='none'
        , ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$y_{Y,0}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('png/Y_residual.png', dpi=300)
    plt.show(block=False)

    plt.figure() 
    plt.plot(yZin, epsilon_expt, color = 'k', marker='o', markerfacecolor='none'
        , ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$y_{Z,0}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('png/Z_residual.png', dpi=300)
    plt.show()

if __name__=="__main__":
    deliverables()