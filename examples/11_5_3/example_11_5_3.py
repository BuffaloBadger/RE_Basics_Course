"""Calculations for Example 11.5.3 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ates
from reb_utils import lls_parameters
from reb_utils import fit_to_SR_data
from reb_utils import Arrhenius_parameters
import os

# global constants available to all functions
# given
V = 3.0 # gal
# known
R = 1.987 # BTU/lbmol/°R

# experimental data as arrays
expt_df = pd.read_csv('example_11_5_3_data.csv')
T = expt_df['T (F)'].to_numpy() + 459.7
Vdot = expt_df['Vdot (gal/min)'].to_numpy()
CA_in = expt_df['CA_in (lbmol/gal)'].to_numpy()
CY_in = expt_df['CY_in (lbmol/gal)'].to_numpy()
CZ_in = expt_df['CZ_in (lbmol/gal)'].to_numpy()
CA_out = expt_df['CA_out (lbmol/gal)'].to_numpy()

# global variables for the current values of the parameters and the 
# index of the current experiment
g_k0_f = float('nan')
g_k0_r = float('nan')
g_E_f = float('nan')
g_E_r = float('nan')
g_iExpt = -1

# CSTR model function
def cstr_model_variables(iExpt):
    # make the index of the current experiment available to the 
    # residuals function
    global g_iExpt
    g_iExpt = iExpt

    # define a guess for the CSTR model variables
    nA_guess = CA_in[iExpt]*Vdot[iExpt]
    nY_guess = CY_in[iExpt]*Vdot[iExpt]
    nZ_guess = CZ_in[iExpt]*Vdot[iExpt]
    initial_guess = np.array([nA_guess, nY_guess, nZ_guess])

    # solve the CSTR design equations
    soln, success, message = solve_ates(cstr_residuals, initial_guess)

    # check for solver issues
    if not success:
        print('')
        print(f'CSTR model function issue for experiment {iExpt}: {message}')
        print('')
        input('Press return to continue of CTRL-C to exit.')
    
    # extract the CSTR model variables
    nA_out = soln[0]
    nY_out = soln[1]
    nZ_out = soln[2]

    # return the CSTR model variables
    return nA_out, nY_out, nZ_out

# CSTR residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA_out = guess[0]
    nY_out = guess[1]
    nZ_out = guess[2]

    # calculate the additional unknowns
    nA_in = CA_in[g_iExpt]*Vdot[g_iExpt]
    nY_in = CY_in[g_iExpt]*Vdot[g_iExpt]
    nZ_in = CZ_in[g_iExpt]*Vdot[g_iExpt]
    CA = nA_out/Vdot[g_iExpt]
    CY = nY_out/Vdot[g_iExpt]
    CZ = nZ_out/Vdot[g_iExpt]
    kf = g_k0_f*np.exp(-g_E_f/R/T[g_iExpt])
    kr = g_k0_r*np.exp(-g_E_r/R/T[g_iExpt])
    r = kf*CA**2 - kr*CY*CZ
    
    # evaluate the residuals
    epsilon_1 = nA_in - nA_out - V*r
    epsilon_2 = nY_in - nY_out + V*r
    epsilon_3 = nZ_in - nZ_out + V*r

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3])

# predicted responses function
def predicted_responses(adj_inputs, beta_f, E_f, beta_r, E_r):
    # make the parameters available to the residuals function
    global g_k0_f, g_E_f, g_k0_r, g_E_r
    g_k0_f = 10**beta_f
    g_E_f = E_f
    g_k0_r = 10**beta_r
    g_E_r = E_r

    # allocate storage for the predicted responses
    CA_out_pred = np.ones_like(CA_out)

    # loop through the experiments
    for iExpt, T_R in enumerate(T):
        # solve the CSTR design equations
        nA_out, nY_out, nZ_out = cstr_model_variables(iExpt)

        # calculate the predicted response
        CA_out_pred[iExpt] = nA_out/Vdot[iExpt]
    
    # return the responses
    return CA_out_pred

# deliverables function
def deliverables():
    # fitting function analysis - define guesses for the parameters
    #par_guess = np.array([0, 1E4, 0, 1E4]) # did not converge
    par_guess = np.array([6.0, 1E4, 6.0, 1E4])

    # combine the adjusted inputs into a matrix
    adj_inputs = np.transpose(np.array([T, Vdot, CA_in, CY_in, CZ_in]))

    # estimate the parameters
    param, param_ci, ff_r_sq, ff_CA_out = fit_to_SR_data(par_guess
            , adj_inputs, CA_out, predicted_responses, use_rel_error=False)

    # extract the results
    ff_k0_f = 10.**param[0]
    ff_k0_f_CI = 10.**param_ci[0,:]
    ff_E_f = param[1]
    ff_E_f_CI = param_ci[1,:]
    ff_k0_r = 10.**param[2]
    ff_k0_r_CI = 10.**param_ci[2,:]
    ff_E_r = param[3]
    ff_E_r_CI = param_ci[3,:]

    # calculate the experiment residuals
    ff_epsilon_expt = CA_out - ff_CA_out

    # linear model analysis - get the data block temperatures
    block_T_values = np.array(expt_df['T (F)'].unique())

    # allocate storage for the results
    kf = np.ones_like(block_T_values)*float('nan')
    kf_CI_lower = np.ones_like(block_T_values)*float('nan')
    kf_CI_upper = np.ones_like(block_T_values)*float('nan')
    kr = np.ones_like(block_T_values)*float('nan')
    kr_CI_lower = np.ones_like(block_T_values)*float('nan')
    kr_CI_upper = np.ones_like(block_T_values)*float('nan')
    model_r_sq = np.ones_like(block_T_values)*float('nan')
    model_plot_data = []

    # loop through the same-temperature data blocks
    for iBlock, blockT in enumerate(block_T_values):
        # extract the same-temperature data block from the full data set
        block_df = expt_df[expt_df['T (F)'] == blockT]

        # extract the data as arrays
        TBlock = block_df['T (F)'].to_numpy() + 459.7
        VdotBlock = block_df['Vdot (gal/min)'].to_numpy()
        CA_inBlock = block_df['CA_in (lbmol/gal)'].to_numpy()
        CY_inBlock = block_df['CY_in (lbmol/gal)'].to_numpy()
        CZ_inBlock = block_df['CZ_in (lbmol/gal)'].to_numpy()
        CA_outBlock = block_df['CA_out (lbmol/gal)'].to_numpy()

        # calculate and save x and y
        tau = V/VdotBlock
        CY_out = CY_inBlock + CA_inBlock - CA_outBlock
        CZ_out = CZ_inBlock + CA_inBlock - CA_outBlock
        y = (CA_inBlock - CA_outBlock)/(tau*CY_out*CZ_out)
        x = CA_outBlock**2/(CY_out*CZ_out)
        model_plot_data.append(x)
        model_plot_data.append(y)

        # fit the linear model to the data
        param, param_ci, r_squared, y_pred = lls_parameters(y, x
                , model_has_intercept=True, use_rel_errors=False)
        model_plot_data.append(y_pred)
        
        # extract and save the results
        kr[iBlock] = -param[0]
        kr_CI_lower[iBlock] = -param_ci[0,0]
        kr_CI_upper[iBlock] = -param_ci[0,1]
        kf[iBlock] = param[1]
        kf_CI_lower[iBlock] = -param_ci[1,0]
        kf_CI_upper[iBlock] = -param_ci[1,1]
        model_r_sq[iBlock] = r_squared

    # fit the Arrhenius expression to the T-kf and T-kr data
    block_T_values = block_T_values + 459.7
    lm_k0_f, lm_k0_f_CI, lm_E_f, lm_E_f_CI, r_sq_Arr_f = Arrhenius_parameters(kf
                , block_T_values, R)
    lm_k0_r, lm_k0_r_CI, lm_E_r, lm_E_r_CI, r_sq_Arr_r = Arrhenius_parameters(kr
                , block_T_values, R)

    # set parameters for calculating the predicted responses
    global g_k0_f, g_E_f, g_k0_r, g_E_r
    g_k0_f = lm_k0_f
    g_E_f = lm_E_f
    g_k0_r = lm_k0_r
    g_E_r = lm_E_r

    # allocate storage for the predicted responses
    lm_CA_out_pred = np.ones_like(CA_out)

    # loop through the experiments
    for iExpt, T_R in enumerate(T):
        # solve the CSTR design equations
        nA_out, nY_out, nZ_out = cstr_model_variables(iExpt)

        # calculate the predicted response
        lm_CA_out_pred[iExpt] = nA_out/Vdot[iExpt]

    # calculate the experiment residuals
    lm_epsilon_expt = CA_out - lm_CA_out_pred

    # calculate the overall r_squared
    CA_mean = np.mean(CA_out)
    ss_res = np.sum(np.square(CA_out - lm_CA_out_pred))
    ss_tot = np.sum(np.square(CA_out - CA_mean))
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
    reb_book_path = '../../../RE_Basics/solutions/ch11_ex3/'
    if not os.path.isdir(reb_book_path):
        # create the folder
        os.makedirs(reb_book_path)

    # tabulate, show and save the fitting results
    fitting_results = [['k0_f', f'{ff_k0_f:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0_f_lower_limit', f'{ff_k0_f_CI[0]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0_f_upper_limit', f'{ff_k0_f_CI[1]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['E_f', f'{ff_E_f:.3g}', 'BTU lbmol^-1^'],
        ['E_f_lower_limit', f'{ff_E_f_CI[0]:.3g}', 'BTU lbmol^-1^'],
        ['E_f_upper_limit', f'{ff_E_f_CI[1]:.3g}', 'BTU lbmol^-1^'],
        ['k0_r', f'{ff_k0_r:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0_r_lower_limit', f'{ff_k0_r_CI[0]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0_r_upper_limit', f'{ff_k0_r_CI[1]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['E_r', f'{ff_E_r:.3g}', 'BTU lbmol^-1^'],
        ['E_r_lower_limit', f'{ff_E_r_CI[0]:.3g}', 'BTU lbmol^-1^'],
        ['E_r_upper_limit', f'{ff_E_r_CI[1]:.3g}', 'BTU lbmol^-1^'],
        ['R_squared', f'{ff_r_sq:.3g}', '']]
    fitting_results_df = pd.DataFrame(fitting_results, columns=['quantity'
            ,'value','units'])
    print(" ")
    print(fitting_results_df)
    print(" ")
    fitting_results_df.to_csv("./csv/example_11_5_3_fit_fcn_results.csv"
            , index=False)
    fitting_results_df.to_csv(reb_book_path + 
        "example_11_5_3_fit_fcn_results.csv", index=False)

    # tabulate, show and save the fitting results
    fitting_results = [['k0_f', f'{lm_k0_f:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0_f_lower_limit', f'{lm_k0_f_CI[0]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0_f_upper_limit', f'{lm_k0_f_CI[1]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['E_f', f'{lm_E_f:.3g}', 'BTU lbmol^-1^'],
        ['E_f_lower_limit', f'{lm_E_f_CI[0]:.3g}', 'BTU lbmol^-1^'],
        ['E_f_upper_limit', f'{lm_E_f_CI[1]:.3g}', 'BTU lbmol^-1^'],
        ['R_squared', f'{r_sq_Arr_f:.3g}', ''],
        ['k0_r', f'{lm_k0_r:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0_r_lower_limit', f'{lm_k0_r_CI[0]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0_r_upper_limit', f'{lm_k0_r_CI[1]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['E_r', f'{lm_E_r:.3g}', 'BTU lbmol^-1^'],
        ['E_r_lower_limit', f'{lm_E_r_CI[0]:.3g}', 'BTU lbmol^-1^'],
        ['E_r_upper_limit', f'{lm_E_r_CI[1]:.3g}', 'BTU lbmol^-1^'],
        ['R_squared', f'{r_sq_Arr_r:.3g}', '']]
    fitting_results_df = pd.DataFrame(fitting_results, columns=['quantity'
            ,'value','units'])
    print(" ")
    print(fitting_results_df)
    print(" ")
    fitting_results_df.to_csv("./csv/example_11_5_3_linear_results.csv"
            , index=False)
    fitting_results_df.to_csv(reb_book_path + 
        "example_11_5_3_linear_results.csv", index=False)

    # create, show, and save Arrhenius plots
    kf_pred = lm_k0_f*np.exp(-lm_E_f/R/block_T_values)
    kr_pred = lm_k0_r*np.exp(-lm_E_r/R/block_T_values)
    plt.figure()
    plt.semilogy(1/block_T_values,kf,color='tab:blue',marker='o', ls='none'
                 , label=f'forward, R$^2$ = {r_sq_Arr_f:.3f}')
    plt.semilogy(1/block_T_values,kf_pred,color='tab:blue')
    plt.semilogy(1/block_T_values,kr,color='tab:orange',marker='o', ls='none'
                 , label=f'reverse, R$^2$ = {r_sq_Arr_r:.3f}')
    plt.semilogy(1/block_T_values,kr_pred,color='tab:orange')
    plt.xlabel('T$^{-1}$ (°R$^{-1}$)')
    plt.ylabel('k (gal lbmol$^{-1}$ min$^{-1}$)')
    plt.xticks(rotation=25)
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_3_Arrhenius_plots.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_3_Arrhenius_plots.png'
            , dpi=300)
    plt.savefig('./pdf/example_11_5_3_Arrhenius_plots.pdf')
    plt.show(block=False)

    # tabulate combined fitting results
    combined_results = []
    combined_results.append(
        ['fitting function', 'k0f', f'{ff_k0_f:.3g} gal/lbmol/min'
            , f'[{ff_k0_f_CI[0]:.3g}, {ff_k0_f_CI[1]:.3g}]'])
    combined_results.append(
        ['', 'Ef', f'{ff_E_f:.0f} BTU/lbmol'
            , f'[{ff_E_f_CI[0]:.0f}, {ff_E_f_CI[1]:.0f}]'])
    combined_results.append(
        ['', 'k0r', f'{ff_k0_r:.3g} gal/lbmol/min'
            , f'[{ff_k0_r_CI[0]:.3g}, {ff_k0_r_CI[1]:.3g}]'])
    combined_results.append(
        ['', 'Er', f'{ff_E_r:.0f} BTU/lbmol'
            , f'[{ff_E_r_CI[0]:.0f}, {ff_E_r_CI[1]:.0f}]'])
    combined_results.append(
        ['','R^2',f'{ff_r_sq:.3f}',''])
    combined_results.append(
        ['linear model', 'k0f', f'{lm_k0_f:.3g} gal/lbmol/min'
            , f'[{lm_k0_f_CI[0]:.3g}, {lm_k0_f_CI[1]:.3g}]'])
    combined_results.append(
        ['', 'Ef', f'{lm_E_f:.0f} BTU/lbmol'
            , f'[{lm_E_f_CI[0]:.0f}, {lm_E_f_CI[1]:.0f}]'])
    combined_results.append(
        ['', 'k0r', f'{lm_k0_r:.3g} gal/lbmol/min'
            , f'[{lm_k0_r_CI[0]:.3g}, {lm_k0_r_CI[1]:.3g}]'])
    combined_results.append(
        ['', 'Er', f'{lm_E_r:.0f} BTU/lbmol'
            , f'[{lm_E_r_CI[0]:.0f}, {lm_E_r_CI[1]:.0f}]'])
    combined_results.append(
        ['','R^2',f'{lm_r_squared:.3f}',''])
    combined_results_df = pd.DataFrame(combined_results
            , columns=['method', 'parameter', 'value', '95% CI'])
    print('')
    print('Combined results')
    print('')
    print(combined_results_df)
    print('')
    combined_results_df.to_csv('./csv/example_11_5_3_combined_results.csv'
            , index=False)
    combined_results_df.to_csv(reb_book_path + 
        'example_11_5_3_combined_results.csv', index=False)

    # generate, show and save a combined parity plot
    plt.figure()
    plt.plot(CA_out, lm_CA_out_pred, color = 'tab:blue', marker='x', ls=''
             , label=f'linear model, R$^2$ = {ff_r_sq:.3f}') 
    plt.plot(CA_out, ff_CA_out, color = 'tab:orange', marker='+', ls=''
             , label=f'fitting function, R$^2$ = {lm_r_squared:.3f}')
    plt.plot([min(CA_out),max(CA_out)],[min(CA_out),max(CA_out)]
             , color = 'k', ls = '-'
             , label = 'Parity Line')
    plt.xlabel('$C_{A,out,meas}$ (lbmol gal$^{-1}$)')
    plt.xlim(left=0)
    plt.ylabel('$C_{A,out,pred}$ (lbmol gal$^{-1}$)')
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_3_parity.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_3_parity.png', dpi=300)
    plt.savefig('./pdf/example_11_5_3_parity.pdf')
    plt.show(block=False)

    # generate show and save residuals plots
    plt.figure() 
    plt.plot(T + 459.7, lm_epsilon_expt, markerfacecolor='none'
             , color = 'tab:blue', marker='x', ls='', label='linear model') 
    plt.plot(T + 459.7, ff_epsilon_expt, markerfacecolor='none'
             , color = 'tab:orange', marker='+', ls='', label='fitting function')
    plt.axhline(y=0, color = 'k')
    plt.xlabel('Temperature (°F)')
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_3_T_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_3_T_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_3_T_residual.pdf')
    plt.show(block=False)

    plt.figure() 
    plt.plot(Vdot, lm_epsilon_expt, markerfacecolor='none', color = 'tab:blue'
             , marker='x', ls='', label='linear model')
    plt.plot(Vdot, ff_epsilon_expt, markerfacecolor='none', color = 'tab:orange'
             , marker='+', ls='', label='fitting function')
    plt.axhline(y=0, color = 'k')
    plt.xlabel('Flow Rate (gal min$^{-1}$)')
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_3_Vdot_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_3_Vdot_residual.png'
            , dpi=300)
    plt.savefig('./pdf/example_11_5_3_Vdot_residual.pdf')
    plt.show(block=False)

    plt.figure()  
    plt.plot(CA_in, lm_epsilon_expt, markerfacecolor='none', color = 'tab:blue'
             , marker='x', ls='', label='linear model')
    plt.plot(CA_in, ff_epsilon_expt, markerfacecolor='none'
             , color = 'tab:orange', marker='+', ls='', label='fitting function')
    plt.axhline(y=0, color = 'k')
    plt.xlabel('Concentration of A (lbmol gal$^{-1}$)')
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_3_CA_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_3_CA_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_3_CA_residual.pdf')
    plt.show(block=False)

    plt.figure()  
    plt.plot(CY_in, lm_epsilon_expt, markerfacecolor='none', color = 'tab:blue'
             , marker='x', ls='', label='linear model')
    plt.plot(CY_in, ff_epsilon_expt, markerfacecolor='none'
             , color = 'tab:orange', marker='+', ls='', label='fitting function')
    plt.axhline(y=0, color = 'k')
    plt.xlabel('Concentration of Y (lbmol gal$^{-1}$)')
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_3_CY_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_3_CY_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_3_CY_residual.pdf')
    plt.show(block=False)

    plt.figure()  
    plt.plot(CZ_in, lm_epsilon_expt, markerfacecolor='none', color = 'tab:blue'
             , marker='x', ls='', label='linear model')
    plt.plot(CZ_in, ff_epsilon_expt, markerfacecolor='none'
             , color = 'tab:orange', marker='+', ls='', label='fitting function')
    plt.axhline(y=0, color = 'k')
    plt.xlabel('Concentration of Z (lbmol gal$^{-1}$)')
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_3_CZ_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_3_CZ_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_3_CZ_residual.pdf')
    plt.show(block=False)

    # combined model plots for linear model
    plot_data = np.transpose(np.array(model_plot_data))
    plt.figure()
    colors = ['tab:blue','tab:orange','tab:green','tab:purple']
    for i, T_R in enumerate(block_T_values):
        ii = 3*i
        x = plot_data[:,ii]
        ym = plot_data[:,ii+1]
        yp = plot_data[:,ii+2]
        plt.plot(x, ym, marker='o', ls='', color=colors[i]
                 , label=f'{T_R - 459.7:.0f} °F, R$^2$ = {model_r_sq[i]:.3f}')
        plt.plot(x, yp, ls='-', color=colors[i])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(title='linear model')
    plt.savefig('./pdf/example_11_5_3_linear_model_plots.pdf')
    plt.savefig('./png/example_11_5_3_linear_model_plots.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_3_linear_model_plots.png', dpi=300)
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    