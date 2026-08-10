"""Calculations for the Class 26 Practice Assignment in REB, The Course"""

# import libraries
import pandas as pd
import numpy as np
from reb_utils import solve_ivodes
from reb_utils import fit_to_SR_data
import matplotlib.pyplot as plt
import os

# global constants available to all functions
# given
D = 2.5 # cm
L = 100.0 # cm
P = 1.0 # atm
K0 = 2.2E23 # atm^3
dH = 200.0 # kJ /mol
# known
Ren = 8.3145E-3 # kJ /mol /K
Rpv = 82.06 # cc atm /mol /K

# experimental data as arrays
expt_df = pd.read_csv('practice_26_data.csv')
VdotAin = expt_df['VFR_A_in'].to_numpy()
VdotYin = expt_df['VFR_Y_in'].to_numpy()
VdotZin = expt_df['VFR_Z_in'].to_numpy()
T = expt_df['T'].to_numpy() + 273.15
alpha = expt_df['Y_to_A'].to_numpy()

# globally variables for the current values of the parameters
global g_k0, g_E, g_T
g_k0 = float('nan')
g_E = float('nan')
g_T = float('nan')

# PFR model function
def pfr_model_variables(VFR_A, VFR_Y, VFR_Z, T, k0, E):
    # make current values of k0, E, T, and tau available to the derivatives function
    global g_k0, g_E, g_T
    g_k0 = k0
    g_E = E
    g_T = T
 
    # define initial values
    ind0 = 0
    dep0 = [P*VFR_A/Rpv/T, P*VFR_Y/Rpv/T, P*VFR_Z/Rpv/T]

    # define the stopping criterion
    f_var = 0
    f_val = L

    # solve the PFR reactor design equations
    z, dep, success, message = solve_ivodes(ind0, dep0, f_var, f_val
            , pfr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not success:
        print('')
        print(f"PFR model solver issue: {message}")
        print('')
        input('press return to continue or CTRL-C to exit')
    
    # extract the individual dependent variables
    nA = dep[0,:]
    nY = dep[1,:]
    nZ = dep[2,:]

    # return the PFR model variables
    return z, nA, nY, nZ

# PFR derivatives function
def pfr_derivatives(ind,dep):
    # extract the dependent variables
    nA = dep[0]
    nY = dep[1]
    nZ = dep[2]

    # calculate the additional unknowns
    k = g_k0*np.exp(-g_E/(Ren*g_T))
    K = K0*np.exp(-dH/Ren/g_T)
    PA = nA*P/(nA + nY + nZ)
    PY = nY*P/(nA + nY + nZ)
    PZ = nZ*P/(nA + nY + nZ)
    r = k*PA*(1 - PY*PZ**3/K/PA)

    # evaluate the derivatives
    dnAdz = -np.pi*D**2/4*r
    dnYdz = np.pi*D**2/4*r
    dnZdz = 3*np.pi*D**2/4*r

    # return the derivatives in an array
    return np.array([dnAdz, dnYdz, dnZdz])

# predicted responses function
def predicted_responses(adj_inputs, beta_guess, E_guess):
    # calculate k from the guess
    k0 = 10**beta_guess

    # allocate storage for the responses
    alpha_pred = np.ones_like(alpha)*float('nan')

    # loop through all of the experiments
    for i, input in enumerate(adj_inputs):
        # solve the PFR design equations
        t, nA, nY, nZ = pfr_model_variables(VdotAin[i], VdotYin[i], VdotZin[i]
            , T[i], k0, E_guess)

        # calculate the response
        alpha_pred[i] = nY[-1]/nA[-1]
    
    # return the predicted responses
    return alpha_pred

# deliverables function
def deliverables():
    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([VdotAin, VdotYin, VdotZin, T]))

    # make a guess for the base-10 log of k0 and for E
    par_guess = [0.0, 40.0]

    # estimate the kinetics parameters
    param, param_ci, r_squared, alpha_pred = fit_to_SR_data(par_guess
        , adj_inputs, alpha, predicted_responses, use_rel_error=True)

    # extract the parameter estimates and their confidence intervals
    beta = param[0]
    beta_CI = param_ci[0,:]
    k0 = 10**beta
    k0_CI = 10**beta_CI 
    E = param[1]
    E_CI = param_ci[1,:]

    # calculate the experiment residuals
    epsilon_expt = alpha - alpha_pred

    # make sure a png folder exists for storing results
    if not os.path.isdir('png'):
        # create the folder
        os.makedirs('./png')
    
    # tabulate, show and save the results to a .csv file
    data = [['k0', f'{k0:.3g}', f'{k0_CI[0]:.3g}'
             , f'{k0_CI[1]:.3g}', 'mol cm^-3^ s^-1^ atm^-1^'],
        ['E', f'{E:.3g}', f'{E_CI[0]:.3g}', f'{E_CI[1]:.3g}', 'kJ mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', float('nan'), float('nan'), ' ']]
    result_df = pd.DataFrame(data, columns=['Parameter','Value','lower_limit'
            , 'upper_limit', 'units'])
    print('')
    print(result_df)
    result_df.to_csv("results.csv", index=False)

    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(alpha, alpha_pred, color = 'b', marker='o'
            , markerfacecolor='none', ls='', label='Data')
    plt.plot([np.min(alpha), np.max(alpha)]
            , [np.min(alpha), np.max(alpha)], color = 'k'
            , label='Parity Line')
    plt.xlabel("Experimental Response")
    plt.ylabel("Model-Predicted Response")
    plt.legend()
    plt.savefig('png/parity.png', dpi=300)
    plt.show(block=False)

    # create, show, and save residuals plots
    plt.figure() 
    plt.plot(VdotAin, epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Inlet A Flow Rate (cm$^3$ s$^{-1}$)")
    plt.ylabel("Residual")
    plt.savefig('png/A_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure() 
    plt.plot(VdotYin, epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Inlet Y Flow Rate (cm$^3$ s$^{-1}$)")
    plt.ylabel("Residual")
    plt.savefig('png/Y_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure() 
    plt.plot(T-273.15, epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Residual")
    plt.savefig('png/T_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure() 
    plt.plot(VdotZin, epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Inlet Z Flow Rate (cm$^3$ s$^{-1}$)")
    plt.ylabel("Residual")
    plt.savefig('png/Z_residuals.png', dpi=300)
    plt.show()

if __name__=="__main__":
    deliverables()
