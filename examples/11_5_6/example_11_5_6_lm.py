"""Calculations for Example 11.5.6 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import lls_parameters
from reb_utils import Arrhenius_parameters
import matplotlib.pyplot as plt
import os

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

# deliverables function
def deliverables():
    # get the data block temperatures
    block_T_values = np.array(expt_df['T'].unique())

    # allocate storage for the results
    lm_k = np.ones_like(block_T_values)*float('nan')
    lm_k_CI_lower = np.ones_like(block_T_values)*float('nan')
    lm_k_CI_upper = np.ones_like(block_T_values)*float('nan')
    lm_model_r_sq = np.ones_like(block_T_values)*float('nan')
    lm_model_plot_data = []

    # loop through the blocks
    for iBlock, blockT in enumerate(block_T_values):
        # extract the same-temperature data block
        block_df = expt_df[expt_df['T'] == blockT]

        # extract the data in the block
        CA0_block = block_df['CA0'].to_numpy()
        t_meas_block = block_df['t_meas'].to_numpy()
        CA_meas_block = block_df['CA_meas'].to_numpy()

        # calculate x and add it to the data for this block
        x = -t_meas_block
        lm_model_plot_data.append(x)

        # calculate y and add it to the data for this block
        y = np.log(CA_meas_block/CA0_block)
        lm_model_plot_data.append(y)

        # fit a straight line through the origin to the x-y data for this block
        param, param_ci, r_squared, y_pred = lls_parameters(y, x
                , model_has_intercept=False, use_rel_errors=False)
        
        # add y_pred to the data for this block
        lm_model_plot_data.append(y_pred)
        
        # save the results for this block
        lm_k[iBlock] = param[0]
        lm_k_CI_lower[iBlock] = param_ci[0][0]
        lm_k_CI_upper[iBlock] = param_ci[0][1]
        lm_model_r_sq[iBlock] = r_squared

    # fit the Arrhenius expression to the T-k data
    T_K = block_T_values + 273.15
    lm_k0, lm_k0_CI, lm_E, lm_E_CI, lm_r_sq_Arr = Arrhenius_parameters(lm_k, T_K
            , R)

    # calculate the predicted responses and residuals
    lm_k_pred = lm_k0*np.exp(-lm_E/R/T_K)
    k = lm_k0*np.exp(-lm_E/R/T)
    lm_CA_pred = CA_0*np.exp(-k*t_meas)
    lm_epsilon_expt = CA_meas - lm_CA_pred

    # calculate the overall r_squared
    CA_mean = np.mean(CA_meas)
    ss_res = np.sum(np.square(CA_meas - lm_CA_pred))
    ss_tot = np.sum(np.square(CA_meas - CA_mean))
    lm_r_sq = 1 - ss_res/ss_tot

    # make sure a results folder and a png folder exist
    if not os.path.isdir('results'):
        # create the folder
        os.makedirs('./results')
    if not os.path.isdir('png'):
        # create the folder
        os.makedirs('./png')

    # tabulate, show, and save the fitting function results
    results = [
        ['k0', lm_k0, lm_k0_CI[0], lm_k0_CI[1], '/min']
        , ['E', lm_E, lm_E_CI[0], lm_E_CI[1], 'kJ/mol']
        , ['R^2 Arr', lm_r_sq_Arr, float('nan'), float('nan'),'']
        , ['R^2 full', lm_r_sq, float('nan'), float('nan'),'']
    ]
    results_df = pd.DataFrame(results)
    results_df.to_csv('results/lm_results.csv', index=False)
    print('')
    print('Arrhenius fitting results:')
    print(results_df)
    print('')

    # tabulate and save the Arrhenius plot data
    results_df = pd.DataFrame({'block_T_values': block_T_values
            , 'lm_k': lm_k, 'lm_k_pred': lm_k_pred})
    results_df.to_csv('results/lm_Arr_plot_data.csv', index=False)
    
    # tabulate, show, and save the model plot results
    results_df = pd.DataFrame({'block_T_values': block_T_values
            , 'lm_k': lm_k, 'lm_k_CI_lower': lm_k_CI_lower
            , 'lm_k_CI_upper': lm_k_CI_upper, 'lm_model_r_sq': lm_model_r_sq})
    results_df.to_csv('results/lm_model_plot_results.csv', index=False)

    # tabulate and save the parity/residual plots data
    results_df = pd.DataFrame({'lm_CA_pred': lm_CA_pred
            , 'lm_epsilon_expt': lm_epsilon_expt})
    results_df.to_csv('results/lm_parity_plot_data.csv', index=False)

    # tabulate and save the model plot data
    plot_data = np.transpose(np.array(lm_model_plot_data))
    results_df = pd.DataFrame(plot_data)
    results_df.to_csv('results/lm_model_plot_data.csv', index=False)

    # plot, show, and save the combined model plots
    plt.figure()
    plt.plot(plot_data[:,0], plot_data[:,1], marker='x', ls='', color='tab:blue'
            , label=f'{T_K[0] - 273.15:.0f} °C')
    plt.plot(plot_data[:,0], plot_data[:,2], ls='-', color='tab:blue')
    plt.plot(plot_data[:,3], plot_data[:,4], marker='x', ls=''
            , color='tab:orange', label=f'{T_K[1] - 273.15:.0f} °C')
    plt.plot(plot_data[:,3], plot_data[:,5], ls='-', color='tab:orange')
    plt.plot(plot_data[:,6], plot_data[:,7], marker='x', ls=''
            , color='tab:green', label=f'{T_K[2] - 273.15:.0f} °C')
    plt.plot(plot_data[:,6], plot_data[:,8], ls='-', color='tab:green')
    plt.plot(plot_data[:,9], plot_data[:,10], marker='x', ls=''
            , color='tab:purple', label=f'{T_K[3] - 273.15:.0f} °C')
    plt.plot(plot_data[:,9], plot_data[:,11], ls='-', color='tab:purple')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.savefig('png/lm_model_plots.png', dpi=300)
    plt.show(block=False)

    # plot, show, and save the Arrhenius plot
    plt.figure()
    plt.semilogy(1/T_K,lm_k,color='tab:blue',marker='o', ls='none')
    plt.semilogy(1/T_K,lm_k_pred,color='k')
    plt.xlabel('T$^{-1}$ (K$^{-1}$)')
    plt.ylabel('k (min$^{-1}$)')
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.savefig('png/lm_Arrhenius_plot.png', dpi=300)
    plt.show(block=False)

    # generate, show, and save the overall parity plot
    plt.figure()
    plt.plot(CA_meas, lm_CA_pred, marker='x', ls='', color='tab:blue'
             ,markerfacecolor='none', label='Data')
    plt.plot([min(CA_meas),max(CA_meas)],[min(CA_meas),max(CA_meas)]
             , color = 'k', ls = '-', label = 'Parity Line')
    plt.xlabel('$C_{A,meas}$ (M)')
    plt.ylabel('$C_{A,pred}$ (M)')
    plt.legend()
    plt.savefig('./png/lm_parity.png', dpi=300)
    plt.show(block=False)

    # generate, show and save the residuals plots
    plt.figure()
    plt.plot(CA_0, lm_epsilon_expt, markerfacecolor='none', marker='x'
             , color='tab:blue', ls='')
    plt.axhline(y=0, color='k')
    plt.xlabel('$C_{A,0}$ (M)')
    plt.ylabel('Residual (M)')
    plt.tight_layout()
    plt.savefig('./png/lm_CA_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure()
    plt.plot(T-273.15, lm_epsilon_expt, markerfacecolor='none', marker='x'
             , color='tab:blue', ls='')
    plt.axhline(y=0, color='k')
    plt.xlabel('T (°C)')
    plt.ylabel('Residual (M)')
    plt.tight_layout()
    plt.savefig('./png/lm_T_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure()
    plt.plot(t_meas, lm_epsilon_expt, markerfacecolor='none', marker='x'
             , color='tab:blue', ls='')
    plt.axhline(y=0, color='k')
    plt.xlabel('$t_{meas}$ (min)')
    plt.ylabel('Residual (M)')
    plt.tight_layout()
    plt.savefig('./png/lm_time_residuals.png', dpi=300)
    plt.show()


# execution command
if __name__ == '__main__':
    deliverables()
    