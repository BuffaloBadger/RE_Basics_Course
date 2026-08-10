"""Calculations for the Class 25 Learning Activity in REB, The Course"""

# import libraries
import pandas as pd
import numpy as np
from reb_utils import solve_ivodes
from reb_utils import fit_to_SR_data
import matplotlib.pyplot as plt
import os

# given and known constants
V = 50.0E-3 # L

# experimental data as arrays
expt_df = pd.read_csv('activity_25_data.csv')
CS0 = expt_df['CS0'].to_numpy()
t_meas = expt_df['t_meas'].to_numpy()
CP_meas = expt_df['CP_meas'].to_numpy()

# globally variables for the current values of the parameters
global g_Vmax, g_Km
g_Vmax = float('nan')
g_Km = float('nan')

# BSTR model function
def bstr_model_variables(CS0, t_meas, Vmax , Km):
    # make current values of Vmax and Km available to the derivatives function
    global g_Vmax, g_Km
    g_Vmax = Vmax
    g_Km = Km
 
    # define initial values
    ind0 = 0
    dep0 = [CS0*V, 0.0]

    # define the stopping criterion
    f_var = 0
    f_val = t_meas

    # solve the BSTR reactor design equations
    t, dep, success, message = solve_ivodes(ind0, dep0, f_var, f_val
            , bstr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model solver issue: {message}")
        print('')
        input('press return to continue or CTRL-C to exit')
    
    # extract the individual dependent variables
    nS = dep[0,:]
    nP = dep[1,:]

    # return the BSTR model variables
    return t, nS, nP

# derivatives function
def bstr_derivatives(ind,dep):
    # extract the dependent variables
    nS = dep[0]

    # calculate the additional unknowns
    CS = nS/V
    r = g_Vmax*CS/(g_Km + CS)

    # evaluate the derivatives
    dnSdt = -r*V
    dnPdt = r*V

    # return the derivatives in an array
    return np.array([dnSdt, dnPdt])

# predicted responses function
def predicted_responses(adj_inputs, Vmax_guess, Km_guess):
    # calculate Vmax and Km from the guesses
    Vmax = 10**Vmax_guess
    Km = 10**Km_guess

    # allocate storage for the responses
    resp = np.ones_like(t_meas)*float('nan')

    # loop through all of the experiments
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs and make T available to all functions
        CS0 = input[0]
        tf = input[1]

        # solve the BSTR design equations
        t, nS, nP = bstr_model_variables(CS0, tf, Vmax, Km)

        # calculate the response
        nPf = nP[-1]
        resp[i] = nPf/V
    
    # return the predicted responses
    return resp

# deliverables function
def deliverables():
    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([CS0, t_meas, CP_meas]))

    # make a guess for the base-10 log of Vmax and Km
    par_guess = [0.0, 0.0]

    # estimate the kinetics parameters
    param, param_ci, r_squared, CP_pred = fit_to_SR_data(par_guess
        , adj_inputs, CP_meas, predicted_responses, use_rel_error=False)

    # extract the parameter estimates and their confidence intervals
    Vmax = 10**param[0]
    Vmax_CI = 10**param_ci[0,:]
    Km = 10**param[1]
    Km_CI = 10**param_ci[1,:]

    # calculate the experiment residuals
    epsilon_expt = CP_meas - CP_pred

    # make sure folders exist for storing results
    if not os.path.isdir('results'):
        # create the folder
        os.makedirs('./results')
    if not os.path.isdir('pdf'):
        # create the folder
        os.makedirs('./pdf')
    if not os.path.isdir('png'):
        # create the folder
        os.makedirs('./png')
    
    # tabulate, show and save the results to a .csv file
    data = [['Vmax', f'{Vmax:.3g}', f'{Vmax_CI[0]:.3g}'
             , f'{Vmax_CI[1]:.3g}', 'mmol L^-1^ min^-1^'],
        ['Km', f'{Km:.3g}', f'{Km_CI[0]:.3g}', f'{Km_CI[1]:.3g}', 'mmol L^-1^'],
        ['R_squared', f'{r_squared:.3g}', float('nan'), float('nan'), ' ']]
    result = pd.DataFrame(data, columns=['Parameter','Value','lower_limit'
            , 'upper_limit', 'units'])
    result.to_csv("results/results.csv", index=False)

    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(CP_meas, CP_pred, color = 'b', marker='o'
            , markerfacecolor='none', ls='', label='Data')
    plt.plot([np.min(CP_meas), np.max(CP_meas)]
            , [np.min(CP_meas), np.max(CP_meas)], color = 'k'
            , label='Parity Line')
    plt.xlabel("Experimental Response (mmol L$^{-1}$)")
    plt.ylabel("Model-Predicted Response (mmol L$^{-1}$)")
    plt.legend()
    plt.savefig('png/parity.png')
    plt.show(block=False)

    # create, show, and save residuals plots
    plt.figure() 
    plt.plot(t_meas, epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Reaction time (min)")
    plt.ylabel("Residual (mmol L$^{-1}$)")
    plt.savefig('png/time_residuals.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(CS0, epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Initial Substrate Concentration (mmol L$^{-1}$)")
    plt.ylabel("Residual (mmol L$^{-1}$)")
    plt.savefig('png/CS0_residuals.png')
    plt.show()

if __name__=="__main__":
    deliverables()
