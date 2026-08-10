"""Calculations for the Class 25 Practice Assignment from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ivodes
from reb_utils import fit_to_SR_data
import os
import matplotlib.pyplot as plt

# global constants available to all functions
# given
V = 2.0 # L
P0 = 1.0 # atm
nZ0 = 0.0 # mol
# known
Ren = 8.314E-3 # kJ/mol/K
Rpv = 0.08206 # L*atm/mol/K

# experimental data as arrays
expt_df = pd.read_csv('practice_25_data.csv')
expt = expt_df['Experiment'].to_numpy()
T = expt_df['T'].to_numpy() + 273.15 # K
gamma = expt_df['AtoB'].to_numpy()
t_meas = expt_df['t'].to_numpy() # min
yZ_meas = expt_df['yZ'].to_numpy()

# global variables for quantities needed by the derivatives function
g_T = float('nan')
g_k0 = float('nan')
g_E = float('nan')

# BSTR model function
def bstr_model_variables(T_K, AtoB, t_meas, k0, E):
    # make the temperature available to the derivatives function
    global g_T, g_k0, g_E
    g_T = T_K
    g_k0 = k0
    g_E = E

    # define the initial values    
    ind_0 = 0
    n0 = P0*V/(Rpv*T_K)
    nB0 = n0/(1.0 + AtoB)
    nA0 = AtoB*nB0
    dep_0 = np.array([nA0, nB0, nZ0])

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
    return t, dep[0,:], dep[1,:], dep[2,:]

# BSTR derivatives function
def bstr_derivatives(ind, dep):
    # extract the individual dependent variables
    nA = dep[0]
    nB = dep[1]
    nZ = dep[2]

    # calculate the additinal unknowns
    k = g_k0 * np.exp( -g_E/(Ren*g_T) )
    CA = nA/V
    CB = nB/V
    r = k * CA *CB

    # evaluate the derivatives
    dnAdt = -r*V
    dnBdt = -r*V
    dnZdt = r*V

    # return the derivatives as an array
    return np.array([dnAdt, dnBdt, dnZdt])

# predicted responses function
def predicted_responses(adj_inputs, beta, E):
    # convert from the base-10 log of k0 to k0
    k0 = 10**beta

    # allocate storage for the predicted responses
    yZ_pred = np.ones_like(yZ_meas)*float('nan')

    # loop through the experiments
    for iExpt, T_K in enumerate(T):
        # solve the BSTR model equations
        t, nA, nB, nZ = bstr_model_variables(T_K, gamma[iExpt]
            , t_meas[iExpt], k0, E)

        # calculate the predicted response
        yZ_pred[iExpt] = nZ[-1]/(nA[-1] + nB[-1] + nZ[-1])

    # return the predicted responses
    return yZ_pred

# deliverables function
def deliverables():
    # fitting function analysis - define a guess for the parameters
    #par_guess = np.array([0.0, 40.0])
    par_guess = np.array([6.0, 40.0])

    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([T, gamma, t_meas]))

    # fit the bstr model to the data
    param, param_ci, r_sq, yZ_pred = fit_to_SR_data(par_guess
            , adj_inputs, yZ_meas, predicted_responses
            , use_rel_error=False)
    
    # extract the results and convert from the base-10 log of k0 to k0
    k0 = 10.**param[0]
    k0_CI = 10.**param_ci[0,:]
    E = param[1]
    E_CI = param_ci[1,:]

    # calculate the experiment residuals
    epsilon_expt = yZ_meas - yZ_pred

    # make sure a results folder and a png folder exist
    if not os.path.isdir('png'):
        # create the folder
        os.makedirs('./png')

    # tabulate, show, and save the fitting function results
    results = [
        ['k0', k0, k0_CI[0], k0_CI[1], '/min']
        , ['E', E, E_CI[0], E_CI[1], 'kJ/mol']
        , ['R-squared', r_sq, '', '', '']
    ]
    results_df = pd.DataFrame(results, columns=['Parameter','Value'
            ,'lower_limit', 'upper_limit', 'units'])
    results_df.to_csv('results.csv', index=False)
    print('')
    print('Fitting function results:')
    print(results_df)
    print('')

    # generate, show, and save the parity plot
    plt.figure()
    plt.plot(yZ_meas, yZ_pred, marker='o', ls='', color='b'
             ,markerfacecolor='none', label='Data')
    plt.plot([min(yZ_meas),max(yZ_meas)],[min(yZ_meas),max(yZ_meas)]
             , color = 'k', ls = '-', label = 'Parity Line')
    plt.xlabel('$y_{Z,meas}$')
    plt.ylabel('$y_{Z,pred}$')
    plt.legend()
    plt.savefig('./png/parity.png', dpi=300)
    plt.show(block=False)

    # generate, show and save the residuals plots
    plt.figure()
    plt.plot(gamma, epsilon_expt, markerfacecolor='none', marker='o'
             , color='b', ls='')
    plt.axhline(y=0, color='k')
    plt.xlabel('Initial A to B Ratio')
    plt.ylabel('Residual')
    plt.tight_layout()
    plt.savefig('./png/A_to_B_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure()
    plt.plot(T-273.15, epsilon_expt, markerfacecolor='none', marker='o'
             , color='b', ls='')
    plt.axhline(y=0, color='k')
    plt.xlabel('T (°C)')
    plt.ylabel('Residual')
    plt.tight_layout()
    plt.savefig('./png/T_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure()
    plt.plot(t_meas, epsilon_expt, markerfacecolor='none', marker='o'
             , color='b', ls='')
    plt.axhline(y=0, color='k')
    plt.xlabel('$t_{meas}$ (min)')
    plt.ylabel('Residual')
    plt.tight_layout()
    plt.savefig('./png/time_residuals.png', dpi=300)
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()