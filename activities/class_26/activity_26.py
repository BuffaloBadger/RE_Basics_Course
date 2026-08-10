"""Calculations for the Class 26 Learning Activity in REB, The Course"""

# import libraries
import pandas as pd
import numpy as np
from reb_utils import solve_ivodes
from reb_utils import fit_to_SR_data
from reb_utils import lls_parameters
from reb_utils import Arrhenius_parameters
import matplotlib.pyplot as plt
import os

# global constants available to all functions
# given
D = 0.3 # dm
L = 10.0 # dm
# known
R = 1.987E-3 # kcal mol^-1^ K^-1^

# experimental data as arrays
expt_df = pd.read_csv('activity_26_data.csv')
T = expt_df['T'].to_numpy() + 273.15
CAin = expt_df['CAin'].to_numpy()
CBin = expt_df['CBin'].to_numpy()
tau = expt_df['tau'].to_numpy()
CA_meas = expt_df['CA_meas'].to_numpy()

# globally variables for the current values of the parameters
global g_k0, g_E, g_T, g_tau
g_k0 = float('nan')
g_E = float('nan')
g_T = float('nan')
g_tau = float('nan')

# PFR model function
def pfr_model_variables(T, CAin, CBin, tau, k0, E):
    # make current values of k0, E, T, and tau available to the derivatives function
    global g_k0, g_E, g_T, g_tau
    g_k0 = k0
    g_E = E
    g_T = T
    g_tau = tau
 
    # define initial values
    ind0 = 0
    Vdot = np.pi*D**2*L/4/tau
    dep0 = [Vdot*CAin, Vdot*CBin, 0.0, 0.0]

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
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]

    # return the PFR model variables
    return z, nA, nB, nY, nZ

# PFR derivatives function
def pfr_derivatives(ind,dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nY = dep[2]
    nZ = dep[3]

    # calculate the additional unknowns
    k = g_k0*np.exp(-g_E/(R*g_T))
    Vdot = np.pi*D**2*L/4/g_tau
    CA = nA/Vdot
    CB = nB/Vdot
    r = k*CA*CB

    # evaluate the derivatives
    dnAdz = -np.pi*D**2/4*r
    dnBdz = -np.pi*D**2/4*r
    dnYdz = np.pi*D**2/4*r
    dnZdz = np.pi*D**2/4*r

    # return the derivatives in an array
    return np.array([dnAdz, dnBdz, dnYdz, dnZdz])

# predicted responses function
def predicted_responses(adj_inputs, beta_guess, E_guess):
    # calculate k from the guess
    k0 = 10**beta_guess

    # allocate storage for the responses
    resp = np.ones_like(CA_meas)*float('nan')

    # loop through all of the experiments
    for i, input in enumerate(adj_inputs):
        # solve the PFR design equations
        z, nA, nB, nY, nZ = pfr_model_variables(T[i], CAin[i], CBin[i], tau[i]
            , k0, E_guess)

        # calculate the response
        Vdot = np.pi*D**2*L/4/tau[i]
        resp[i] = nA[-1]/Vdot
    
    # return the predicted responses
    return resp

# deliverables function
def deliverables():
    # fitting function analysis - combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([T, CAin, CBin, tau]))

    # make a guess for the base-10 log of k0 and for E
    #par_guess = [0.0, 10.0]
    par_guess = [6.0, 10.0]

    # estimate the kinetics parameters
    param, param_ci, r_squared, CA_pred = fit_to_SR_data(par_guess
        , adj_inputs, CA_meas, predicted_responses, use_rel_error=False)

    # extract the parameter estimates and their confidence intervals
    beta = param[0]
    beta_CI = param_ci[0,:]
    k0 = 10**beta
    k0_CI = 10**beta_CI 
    E = param[1]
    E_CI = param_ci[1,:]

    # calculate the experiment residuals
    epsilon_expt = CA_meas - CA_pred

    # make sure folders exist for storing results
    if not os.path.isdir('results'):
        # create the folder
        os.makedirs('./results')
    if not os.path.isdir('png'):
        # create the folder
        os.makedirs('./png')
    
    # tabulate, show and save the results to a .csv file
    data = [['k0', f'{k0:.3g}', f'{k0_CI[0]:.3g}'
             , f'{k0_CI[1]:.3g}', 'L mol^-1^ min^-1^'],
        ['E', f'{E:.3g}', f'{E_CI[0]:.3g}', f'{E_CI[1]:.3g}', 'cal mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', float('nan'), float('nan'), ' ']]
    result_df = pd.DataFrame(data, columns=['Parameter','Value','lower_limit'
            , 'upper_limit', 'units'])
    print('')
    print('Fitting function results:')
    print(result_df)
    result_df.to_csv("results/fitting_function_results.csv", index=False)

    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(CA_meas, CA_pred, color = 'b', marker='o'
            , markerfacecolor='none', ls='', label='Data')
    plt.plot([np.min(CA_meas), np.max(CA_meas)]
            , [np.min(CA_meas), np.max(CA_meas)], color = 'k'
            , label='Parity Line')
    plt.xlabel("Experimental Response (mol L$^{-1}$)")
    plt.ylabel("Model-Predicted Response (mol L$^{-1}$)")
    plt.legend()
    plt.savefig('png/fitting_function_parity.png')
    plt.show(block=False)

    # create, show, and save residuals plots
    plt.figure() 
    plt.plot(CAin, epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Feed CA (mol L$^{-1}$)")
    plt.ylabel("Residual (mol L$^{-1}$)")
    plt.savefig('png/ff_CA_residuals.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(CBin, epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Feed CB (mol L$^{-1}$)")
    plt.ylabel("Residual (mol L$^{-1}$)")
    plt.savefig('png/ff_CB_residuals.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(T-273.15, epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Residual (mol L$^{-1}$)")
    plt.savefig('png/ff_T_residuals.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(tau, epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Space Time (min)")
    plt.ylabel("Residual (mmol L$^{-1}$)")
    plt.savefig('png/ff_tau_residuals.png')
    plt.show(block=False)

    # linear model analysis - get block temperatures
    block_T_values = np.array(expt_df['T'].unique()) + 273.15

    # allocate storage for the results
    k = np.ones_like(block_T_values)*float('nan')
    k_CI_lower = np.ones_like(block_T_values)*float('nan')
    k_CI_upper = np.ones_like(block_T_values)*float('nan')
    r_sq = np.ones_like(block_T_values)*float('nan')

    # start the model plot
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple']
    plt.figure()

    # loop through the blocks
    for iBlock, blockT in enumerate(block_T_values):
        # extract the same-temperature data block
        block_df = expt_df[expt_df['T'] == blockT - 273.15]

        # extract the data in the block
        CAin_block = block_df['CAin'].to_numpy()
        CBin_block = block_df['CBin'].to_numpy()
        tau_block = block_df['tau'].to_numpy()
        CA_meas_block = block_df['CA_meas'].to_numpy()

        # calculate x, y
        x = -tau_block
        y = 1/(CAin_block - CBin_block) * np.log(CAin_block*(CBin_block 
                - CAin_block + CA_meas_block)/(CBin_block*CA_meas_block))

        # fit a straight line through the origin to the x-y data for this block
        param, param_ci, r_squared, y_pred = lls_parameters(y, x
                , model_has_intercept=False, use_rel_errors=False)
        
        # add the data to the model plot
        plt.plot(x, y, color = colors[iBlock], marker='o'
                , markerfacecolor='none', ls=''
                , label=f'T = {blockT-273.15:.1f} °C')
        plt.plot(x, y_pred, colors[iBlock], ls='-')
        
        # save the results for this block
        k[iBlock] = param[0]
        k_CI_lower[iBlock] = param_ci[0][0]
        k_CI_upper[iBlock] = param_ci[0][1]
        r_sq[iBlock] = r_squared
    
    # finish the model plot
    plt.xlabel("x (min)")
    plt.ylabel("y (L mmol$^{-1}$)")
    plt.legend()
    plt.savefig('png/model_plots.png')
    plt.show(block=False)

    # save the model plot fitting data
    results_df = pd.DataFrame({'T': block_T_values - 273.15, 'k': k
        , 'k_CI_lower': k_CI_lower, 'k_CI_upper': k_CI_upper
        , 'R_squared': r_sq})
    print('')
    print('Model plot fitting results:')
    print(results_df)
    results_df.to_csv("results/model_plot_results.csv", index=False)

    # fit the Arrhenius expression to the T-k data
    lm_k0, lm_k0_CI, lm_E, lm_E_CI, lm_r_sq_Arr = Arrhenius_parameters(k
        , block_T_values, R)

    
    # create, show, and save the Arrhenius plot
    y_pred = lm_k0*np.exp(-lm_E/R/block_T_values)
    plt.figure()
    plt.semilogy(1/block_T_values,k,color='b', markerfacecolor='none'
                , marker='o', ls='none')
    plt.semilogy(1/block_T_values,y_pred,color='k')
    plt.xlabel('T$^{-1}$ (K$^{-1}$)')
    plt.ylabel('k (L mol$^{-1}$ min$^{-1}$)')
    plt.xticks(rotation=25)
    plt.legend(title=f'R$^2$ = {lm_r_sq_Arr:.3f}')
    plt.tight_layout()
    plt.savefig('png/arrhenius_plot.png')
    plt.show(block=False)

    # calculate the predicted responses and residuals
    lm_CA_pred = predicted_responses(adj_inputs, np.log10(lm_k0), lm_E)
    lm_epsilon_expt = CA_meas - lm_CA_pred

    # calculate the overall r_squared
    CA_mean = np.mean(CA_meas)
    ss_res = np.sum(np.square(CA_meas - lm_CA_pred))
    ss_tot = np.sum(np.square(CA_meas - CA_mean))
    lm_r_sq = 1 - ss_res/ss_tot
    
    # tabulate, show and save the results to a .csv file
    data = [['lm_k0', f'{lm_k0:.3g}', f'{lm_k0_CI[0]:.3g}'
             , f'{lm_k0_CI[1]:.3g}', 'L mol^-1^ min^-1^'],
        ['E', f'{lm_E:.3g}', f'{lm_E_CI[0]:.3g}', f'{lm_E_CI[1]:.3g}', 'kJ mol^-1^'],
        ['R_squared', f'{lm_r_sq:.3g}', float('nan'), float('nan'), ' ']]
    result_df = pd.DataFrame(data, columns=['Parameter','Value','lower_limit'
            , 'upper_limit', 'units'])
    print('')
    print('Arrhenius expression fitting results:')
    print(result_df)
    result_df.to_csv("results/linear_model_results.csv", index=False)

    # create, show, and save parity and residuals plots for the linear model
    plt.figure()
    plt.plot(CA_meas, lm_CA_pred, color = 'b', marker='o'
            , markerfacecolor='none', ls='', label='Data')
    plt.plot([np.min(CA_meas), np.max(CA_meas)]
            , [np.min(CA_meas), np.max(CA_meas)], color = 'k'
            , label='Parity Line')
    plt.xlabel("Experimental Response (mol L$^{-1}$)")
    plt.ylabel("Model-Predicted Response (mol L$^{-1}$)")
    plt.legend()
    plt.savefig('png/linear_model_parity.png')
    plt.show(block=False)

    plt.figure()
    plt.plot(CAin, lm_epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Feed CA (mol L$^{-1}$)")
    plt.ylabel("Residual (mol L$^{-1}$)")
    plt.savefig('png/lm_CA_residuals.png')
    plt.show(block=False)

    plt.figure()
    plt.plot(CBin, lm_epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Feed CB (mol L$^{-1}$)")
    plt.ylabel("Residual (mol L$^{-1}$)")
    plt.savefig('png/lm_CB_residuals.png')
    plt.show(block=False)

    plt.figure()
    plt.plot(T-273.15, lm_epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Residual (mol L$^{-1}$)")
    plt.savefig('png/lm_T_residuals.png')
    plt.show(block=False)

    plt.figure()
    plt.plot(tau, lm_epsilon_expt, color = 'b', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Space Time (min)")
    plt.ylabel("Residual (mmol L$^{-1}$)")
    plt.savefig('png/lm_tau_residuals.png')
    plt.show()

if __name__=="__main__":
    deliverables()
