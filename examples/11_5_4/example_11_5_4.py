"""Calculations for Example 11.5.4 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes
from reb_utils import lls_parameters
from reb_utils import fit_to_SR_data
from reb_utils import Arrhenius_parameters
import os

# global constants available to all functions
# given
V = 500.0 # cc
# known
Re = 1.987E-3 # kcal/mol/K
Rpv = 82.06 # cc atm/mol/K

# experimental data as arrays
expt_df = pd.read_csv('example_11_5_4_data.csv')
T = expt_df['T'].to_numpy() + 273.15
PA0 = expt_df['PA0'].to_numpy()
PB0 = expt_df['PB0'].to_numpy()
t_meas = expt_df['t_meas'].to_numpy()
fA = expt_df['fA'].to_numpy()

# global variables for the current values of k0, E, and T
g_k0 = float('nan')
g_E = float('nan')
g_T = float('nan')

# BSTR model function
def bstr_model_variables(T, PA0, PB0, tf):
    # make T available to the derivatives function
    global g_T
    g_T = T

    # set initial values and stopping criterion
    t0 = 0
    nA0 = PA0*V/Rpv/T
    nB0 = PB0*V/Rpv/T
    dep0 = np.array([nA0, nB0, 0.0, 0.0])
    stop_var = 0

    # solve the BSTR reactor design equations
    t, dep, success, message = solve_ivodes(t0, dep0, stop_var, tf
            , bstr_derivatives,True)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")
    
    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]

    # return the profiles
    return t, nA, nB, nY, nZ

# BSTR derivatives function
def bstr_derivatives(t,dep):
    # get the dependent variables that are needed
    nA = dep[0]
    nB = dep[1]

    # calculate the partial pressures
    PA = nA*Rpv*g_T/V
    PB = nB*Rpv*g_T/V

    # calculate the rate
    k = g_k0*np.exp(-g_E/Re/g_T)
    r = k*PA*PB

    # calculate the time derivatives of the dependent variables
    dnAdt = -r*V
    dnBdt = -r*V
    dnYdt = r*V
    dnZdt = r*V

    # return an array containing the derivatives
    ddt = np.array([dnAdt, dnBdt, dnYdt, dnZdt])
    return ddt

# predicted responses function
def predicted_responses(adj_inputs, k0, E):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    resp = np.zeros(nExpts)

    # make the rate coefficient available to other functions
    global g_k0, g_E
    g_k0 = k0
    g_E = E

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # get the other experimental inputs
        T = input[0]
        PA0 = input[1]
        PB0 = input[2]
        tf = input[3]

        # solve the reactor design equations
        t, nA, nB, nY, nZ = bstr_model_variables(T, PA0, PB0, tf)
        
        # calculate the model-predicted response
        nAf = nA[-1]
        nA0 = PA0*V/Rpv/T
        fA = (nA0 - nAf)/nA0
        resp[i] = fA

    # return the responses
    return resp

# deliverables function
def deliverables():
    # fitting function analysis - guess the parameters
    par_guess = [0.0, 10.0]

    # combine the adjusted inputs into a matrix
    adj_inputs = np.transpose(np.array([T, PA0, PB0, t_meas]))

    # estimate the parameters
    param, param_ci, ff_r_squared, ff_fA_pred = fit_to_SR_data(par_guess
            , adj_inputs, fA, predicted_responses, use_rel_error=False)
    
    # extract the results
    ff_k0 = param[0]
    ff_k0_CI = param_ci[0,:]
    ff_E = param[1]
    ff_E_CI = param_ci[1,:]

    # calculate the experiment residuals
    ff_epsilon_expt = fA - ff_fA_pred

    # analysis using a linearized model - get the data block temperatures
    block_T_values = np.array(expt_df['T'].unique())

    # allocate storage for the results
    lm_k = np.ones_like(block_T_values)*float('nan')
    lm_k_CI_lower = np.ones_like(block_T_values)*float('nan')
    lm_k_CI_upper = np.ones_like(block_T_values)*float('nan')
    lm_r_sq = np.ones_like(block_T_values)*float('nan')
    lm_model_plot_data = []

    # loop through the same-temperature data blocks
    for iBlock, blockT in enumerate(block_T_values):
        # extract the same-temperature data block from the full data set
        block_df = expt_df[expt_df['T'] == blockT]

        # extract the block data as arrays
        PA0_Block = block_df['PA0'].to_numpy()
        PB0_Block = block_df['PB0'].to_numpy()
        t_meas_Block = block_df['t_meas'].to_numpy()
        fA_Block = block_df['fA'].to_numpy()

        # allocate storage for x and y
        x = np.ones_like(PA0_Block)*float('nan')
        y = np.ones_like(PA0_Block)*float('nan')

        # calculate x and y
        for iExpt, tf in enumerate(t_meas_Block):
            nA0 = PA0_Block[iExpt]*V/Rpv/(blockT + 273.15)
            nB0 = PB0_Block[iExpt]*V/Rpv/(blockT + 273.15)
            nA = nA0*(1-fA_Block[iExpt])
            x[iExpt] = -tf*(Rpv*(blockT + 273.15))**2/V
            if nA0 == nB0:
                y[iExpt] = 1/nA0 - 1/nA
            else:
                y[iExpt] = 1/(nA0 - nB0)*np.log(nA0*(nB0 - nA0 + nA)/(nB0*nA))

        # fit the linear model to the data
        param, param_ci, r_squared, y_pred = lls_parameters(y, x
                , model_has_intercept=False, use_rel_errors=False)
        
        lm_model_plot_data.append(x)
        lm_model_plot_data.append(y)
        lm_model_plot_data.append(y_pred)
        
        # extract the rate coefficient
        lm_k[iBlock] = param[0]
        lm_k_CI_lower[iBlock] = param_ci[0,0]
        lm_k_CI_upper[iBlock] = param_ci[0,1]
        lm_r_sq[iBlock] = r_squared

    # fit the Arrhenius expression to the T-k data
    block_T_values = block_T_values + 273.15
    lm_k0, lm_k0_CI, lm_E, lm_E_CI, lm_r_sq_Arr = Arrhenius_parameters(lm_k
                , block_T_values, Re)

    # allocate storage for the predicted responses
    lm_fA_pred = np.ones_like(fA)

    # set parameters for calculating the predicted responses
    global g_k0, g_E
    g_k0 = lm_k0
    g_E = lm_E

    # loop through the experiments
    for iExpt, T_K in enumerate(T):
        # solve the CSTR design equations
        t, nA, nB, nY, nZ = bstr_model_variables(T_K, PA0[iExpt], PB0[iExpt]
                ,t_meas[iExpt])

        # calculate the predicted response
        nA0 = PA0[iExpt]*V/Rpv/T_K
        lm_fA_pred[iExpt] = (nA0 - nA[-1])/nA0

    # calculate the experiment residuals
    lm_epsilon_expt = fA - lm_fA_pred

    # calculate the overall r_squared
    fA_mean = np.mean(fA)
    ss_res = np.sum(np.square(fA - lm_fA_pred))
    ss_tot = np.sum(np.square(fA - fA_mean))
    lm_r_squared = 1 - ss_res/ss_tot

    # make sure folders exist for storing results
    if not os.path.isdir('csv'):
        # create the folder
        os.makedirs('./csv')
    if not os.path.isdir('pdf'):
        # create the folder
        os.makedirs('./pdf')
    if not os.path.isdir('png'):
        # create the folder
        os.makedirs('./png')
    reb_book_path = '../../../RE_Basics/solutions/ch11_ex4/'
    if not os.path.isdir(reb_book_path):
        # create the folder
        os.makedirs(reb_book_path)

    # tabulate, show and save the fitting function results
    fitting_results = [['k0', f'{ff_k0:.3g}', 'mol cm^-3^ min^-1^ atm^-2^'],
        ['k0_lower_limit', f'{ff_k0_CI[0]:.3g}', 'mol cm^-3^ min^-1^ atm^-2^'],
        ['k0_upper_limit', f'{ff_k0_CI[1]:.3g}', 'mol cm^-3^ min^-1^ atm^-2^'],
        ['E', f'{ff_E:.3g}', 'kcal mol^-1^'],
        ['E_lower_limit', f'{ff_E_CI[0]:.3g}', 'kcal mol^-1^'],
        ['E_upper_limit', f'{ff_E_CI[1]:.3g}', 'kcal mol^-1^'],
        ['R-squared', f'{ff_r_squared:.3g}', '']]
    results_df = pd.DataFrame( fitting_results
            , columns=['quantity','value','units'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv("./csv/example_11_5_4_fit_fcn_results.csv", index=False)
    results_df.to_csv(reb_book_path + "example_11_5_4_fit_fcn_results.csv"
        , index=False)

    # tabulate, show and save the linear model results
    fitting_results = [['k0', f'{lm_k0:.3g}', 'mol cm^-3^ min^-1^ atm^-2^'],
        ['k0_lower_limit', f'{lm_k0_CI[0]:.3g}', 'mol cm^-3^ min^-1^ atm^-2^'],
        ['k0_upper_limit', f'{lm_k0_CI[1]:.3g}', 'mol cm^-3^ min^-1^ atm^-2^'],
        ['E', f'{lm_E:.3g}', 'kcal mol^-1^'],
        ['E_lower_limit', f'{lm_E_CI[0]:.3g}', 'kcal mol^-1^'],
        ['E_upper_limit', f'{lm_E_CI[1]:.3g}', 'kcal mol^-1^'],
        ['R_squared', f'{lm_r_squared:.3g}', '']]
    fitting_results_df = pd.DataFrame(fitting_results, columns=['quantity'
            ,'value','units'])
    print(" ")
    print(fitting_results_df)
    print(" ")
    fitting_results_df.to_csv("./csv/example_11_5_4_Arrhenius_results.csv"
        , index=False)
    fitting_results_df.to_csv(reb_book_path + 
        "example_11_5_4_Arrhenius_results.csv", index=False)
        
    # create, show, and save parity plots
    plt.figure() 
    plt.plot(fA, lm_fA_pred, color = 'tab:blue', marker='x', ls=''
             , label=f'linear model, R$^2$ = {lm_r_squared:.3f}')
    plt.plot(fA, ff_fA_pred, color='tab:orange', marker='+', ls=''
             , label=f'fitting function, R$^2$ = {ff_r_squared:.3f}')
    plt.plot([min(fA),max(fA)],[min(fA),max(fA)],
        color = 'k', ls = '-', label='Parity Line')
    plt.xlabel("Measured Conversion")
    plt.ylabel("Predicted Conversion")
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_4_parity.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_4_parity.png', dpi=300)
    plt.savefig('./pdf/example_11_5_4_parity.pdf')
    plt.show(block=False)

    # create, show and save residuals plots
    plt.figure()
    plt.plot(PA0, lm_epsilon_expt, markerfacecolor='none', color='tab:blue'
             , marker='x', ls='', label='linear model')
    plt.plot(PA0, ff_epsilon_expt, markerfacecolor='none', color = 'tab:orange'
             , marker='+', ls='', label='fitting function')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$P_{A,0}$ (atm)")
    plt.ylabel("Residual")
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_4_PA0_residuals.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_4_PA0_residuals.png', dpi=300)
    plt.savefig('./pdf/example_11_5_4_PA0_residuals.pdf')
    plt.show(block=False)

    plt.figure()
    plt.plot(PB0, lm_epsilon_expt, markerfacecolor='none', color='tab:blue'
             , marker='x', ls='', label='linear model')
    plt.plot(PB0, ff_epsilon_expt, markerfacecolor='none', color = 'tab:orange'
             , marker='+', ls='', label='fitting function')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$P_{B,0}$ (atm)")
    plt.ylabel("Residual")
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_4_PB0_residuals.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_4_PB0_residuals.png', dpi=300)
    plt.savefig('./pdf/example_11_5_4_PB0_residuals.pdf')
    plt.show(block=False)

    plt.figure()
    plt.plot(t_meas, lm_epsilon_expt, markerfacecolor='none', color='tab:blue'
             , marker='x', ls='', label='linear model')
    plt.plot(t_meas, ff_epsilon_expt, markerfacecolor='none', ls=''
             , color = 'tab:orange', marker='+', label='fitting function')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$t_{meas}$ (min)")
    plt.ylabel("Residual")
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_4_t_meas_residuals.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_4_t_meas_residuals.png', dpi=300)
    plt.savefig('./pdf/example_11_5_4_t_meas_residuals.pdf')
    plt.show(block=False)

    plt.figure() 
    plt.plot(T-273.15, lm_epsilon_expt, markerfacecolor='none', color='tab:blue'
             , marker='x', ls='', label='linear model')
    plt.plot(T-273.15, ff_epsilon_expt, markerfacecolor='none', ls=''
             , color = 'tab:orange', marker='+', label='fitting function')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("T (°C)")
    plt.ylabel("Residual")
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_4_T_residuals.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_4_T_residuals.png', dpi=300)
    plt.savefig('./pdf/example_11_5_4_T_residuals.pdf')
    plt.show(block=False)

    # create, show, and save model plots
    plot_data = np.transpose(np.array(lm_model_plot_data))
    colors = ['tab:blue','tab:orange','tab:green','tab:purple']
    plt.figure()
    for i, T_C in enumerate(block_T_values):
        ii = 3*i
        x = plot_data[:,ii]
        ym = plot_data[:,ii+1]
        yp = plot_data[:,ii+2]
        plt.plot(x, ym, marker='o', ls='', color=colors[i]
                 , label=f'{T_C:.0f} °C')
        plt.plot(x, yp, ls='-', color=colors[i])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(title='linear model')
    plt.savefig('./pdf/example_11_5_4_linear_model_plots.pdf')
    plt.savefig('./png/example_11_5_4_linear_model_plots.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_4_linear_model_plots.png'
        , dpi=300)
    plt.show(block=False)    

    # create, show, and save Arrhenius plot
    k_pred = lm_k0*np.exp(-lm_E/Re/block_T_values)
    plt.figure()
    plt.semilogy(1/block_T_values,lm_k,color='k',marker='o', ls='none')
    plt.semilogy(1/block_T_values,k_pred,color='k')
    plt.xlabel('T$^{-1}$ (°R$^{-1}$)')
    plt.ylabel('k (gal lbmol$^{-1}$ min$^{-1}$)')
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.savefig('./png/example_11_5_4_Arrhenius_plot.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_4_Arrhenius_plot.png', dpi=300)
    plt.savefig('./pdf/example_11_5_4_Arrhenius_plot.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    