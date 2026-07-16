"""Fitting Function Calculations for Example 11.5.6 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ivodes
from reb_utils import fit_to_SR_data
import os
import matplotlib.pyplot as plt


# global constants available to all functions
# given
V = 1.0 # L
# known
R = 8.314E-3 # kJ/mol/K

# experimental data as arrays
expt_df = pd.read_csv('example_11_5_6_data.csv')
expt = expt_df['Experiment'].to_numpy()
T = expt_df['T'].to_numpy() + 273.15 # K
CA_0 = expt_df['CA0'].to_numpy() # mol /L
t_meas = expt_df['t_meas'].to_numpy() # min
CA_meas = expt_df['CA_meas'].to_numpy() # mol /L

# global variables for quantities needed by the derivatives function
g_T = float('nan')
g_k0 = float('nan')
g_E = float('nan')

# BSTR model function
def bstr_model_variables(T, k0, E, CA0, t_meas):
    # make the temperature available to the derivatives function
    global g_T, g_k0, g_E
    g_T = T
    g_k0 = k0
    g_E = E

    # define the initial values
    ind_0 = 0
    dep_0 = np.array([CA0*V, 0.0])

    # define the stopping criterion
    f_var = 0
    f_val = t_meas

    # solve the BSTR model equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            , bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f'BSTR model issue: {message}')
        print('')
        input('Press return to continue of CTRL-C to exit.')
    
    # return the BSTR model variables
    return t, dep[0,:], dep[1,:]

# BSTR derivatives function
def bstr_derivatives(ind, dep):
    # extract the individual dependent variables
    nA = dep[0]
    nZ = dep[1]

    # calculate the additinal unknowns
    k = g_k0 * np.exp( -g_E/(R*g_T) )
    CA = nA/V
    r = k * CA

    # evaluate and return the derivatives
    dnAdt = -r*V
    dnBdt = r*V
    return np.array([dnAdt, dnBdt])

# predicted responses function
def predicted_responses(adj_inputs, beta, E):
    # convert from the base-10 log of k0 to k0
    k0 = 10**beta

    # allocate storage for the predicted responses
    CA_pred = np.ones_like(T)*float('nan')

    # loop through the experiments
    for iExpt, T_K in enumerate(T):
        # solve the BSTR model equations
        t, nA, nZ = bstr_model_variables(T_K, k0, E, CA_0[iExpt], t_meas[iExpt])

        # calculate the predicted response
        CA_pred[iExpt] = nA[-1]/V

    # return the predicted responses
    return CA_pred

# deliverables function
def deliverables():
    # fitting function analysis - define a guess for the parameters
    #par_guess = np.array([0.0, 40.0])
    par_guess = np.array([4.0, 40.0])

    # combine the adjusted inputs as a matrix
    adjusted_inputs = np.transpose(np.array([T, CA_0, t_meas]))

    # fit the bstr model to the data
    param, param_ci, ff_r_sq, ff_CA_pred = fit_to_SR_data(par_guess
            , adjusted_inputs, CA_meas, predicted_responses
            , use_rel_error=False)
    
    # extract the results and convert from the base-10 log of k0 to k0
    ff_k0 = 10.**param[0]
    ff_k0_CI = 10.**param_ci[0,:]
    ff_E = param[1]
    ff_E_CI = param_ci[1,:]

    # calculate the experiment residuals
    ff_epsilon_expt = CA_meas - ff_CA_pred

    # make sure a results folder and a png folder exist
    if not os.path.isdir('results'):
        # create the folder
        os.makedirs('./results')
    if not os.path.isdir('png'):
        # create the folder
        os.makedirs('./png')

    # tabulate, show, and save the fitting function results
    results = [
        ['k0', ff_k0, ff_k0_CI[0], ff_k0_CI[1], '/min']
        , ['E', ff_E, ff_E_CI[0], ff_E_CI[1], 'kJ/mol']
        , ['R-squared', ff_r_sq, float('nan'), float('nan'),'']
    ]
    results_df = pd.DataFrame(results)
    results_df.to_csv('results/ff_results.csv', index=False)
    print('')
    print('Fitting function results:')
    print(results_df)
    print('')

    # tabulate and save the fitting function plot data
    results_df = pd.DataFrame({'ff_CA_pred': ff_CA_pred
            , 'ff_epsilon_expt': ff_epsilon_expt})
    results_df.to_csv('results/ff_plot_data.csv', index=False)

    # generate, show, and save the parity plot
    plt.figure()
    plt.plot(CA_meas, ff_CA_pred, marker='x', ls='', color='tab:blue'
             ,markerfacecolor='none', label='Data')
    plt.plot([min(CA_meas),max(CA_meas)],[min(CA_meas),max(CA_meas)]
             , color = 'k', ls = '-', label = 'Parity Line')
    plt.xlabel('$C_{A,meas}$ (M)')
    plt.ylabel('$C_{A,pred}$ (M)')
    plt.legend()
    plt.savefig('./png/ff_parity.png', dpi=300)
    plt.show(block=False)

    # generate, show and save the residuals plots
    plt.figure()
    plt.plot(CA_0, ff_epsilon_expt, markerfacecolor='none', marker='x'
             , color='tab:blue', ls='')
    plt.axhline(y=0, color='k')
    plt.xlabel('$C_{A,0}$ (M)')
    plt.ylabel('Residual (M)')
    plt.tight_layout()
    plt.savefig('./png/ff_CA_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure()
    plt.plot(T-273.15, ff_epsilon_expt, markerfacecolor='none', marker='x'
             , color='tab:blue', ls='')
    plt.axhline(y=0, color='k')
    plt.xlabel('T (°C)')
    plt.ylabel('Residual (M)')
    plt.tight_layout()
    plt.savefig('./png/ff_T_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure()
    plt.plot(t_meas, ff_epsilon_expt, markerfacecolor='none', marker='x'
             , color='tab:blue', ls='')
    plt.axhline(y=0, color='k')
    plt.xlabel('$t_{meas}$ (min)')
    plt.ylabel('Residual (M)')
    plt.tight_layout()
    plt.savefig('./png/ff_time_residuals.png', dpi=300)
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    