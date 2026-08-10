"""Calculations for the Class 29 Learning Activity from REB, The Course"""

# import libraries
import pandas as pd
import numpy as np
from reb_utils import lls_parameters
from reb_utils import Arrhenius_parameters
import matplotlib.pyplot as plt
import os

# global constants available to all functions
# given
V = 100.0 # cm^3
P0 = 6.0 # atm
# known
Rpv = 82.06 # cm^3 atm/mol/K
Ren = 1.987E-3 # kcal/mol

# experimental data as arrays
expt_df = pd.read_csv('activity_29_data.csv')
T = expt_df['T'].to_numpy() + 273.15
PA0 = expt_df['PA0'].to_numpy()
tf = expt_df['t_meas'].to_numpy()
Pf = expt_df['P_meas'].to_numpy()

# deliverables function
def deliverables():
    # get the data block temperatures
    block_T_values = np.array(expt_df['T'].unique()) + 273.15

    # allocate storage for the results
    k = np.ones_like(block_T_values)*float('nan')
    k_CI_lower = np.ones_like(block_T_values)*float('nan')
    k_CI_upper = np.ones_like(block_T_values)*float('nan')
    model_r_sq = np.ones_like(block_T_values)*float('nan')

    # make sure a png folder exists
    if not os.path.isdir('png'):
        os.makedirs('./png')

    # set up the model plots
    colors = ['tab:blue', 'tab:orange', 'tab:green']
    markers = ['o', 'x', '+']
    plt.figure()

    # loop through the blocks
    for iBlock, blockT in enumerate(block_T_values):
        # extract the same-temperature data block
        block_df = expt_df[expt_df['T'] + 273.15 == blockT]

        # extract the data in the block
        PA0_block = block_df['PA0'].to_numpy()
        t_meas_block = block_df['t_meas'].to_numpy()
        P_meas_block = block_df['P_meas'].to_numpy()

        # calculate molar amounts and partial pressures
        nA0 = PA0_block*V/(Rpv*blockT)
        nB0 = (P0 - PA0_block)*V/(Rpv*blockT)
        xi = nA0 + nB0 - P_meas_block*V/(Rpv*blockT)
        nA = nA0 - xi
        nB = nB0 - xi
        PA = nA*Rpv*blockT/V
        PB = nB*Rpv*blockT/V

        # replace any negative partial pressures with zero
        PA[PA<0] = 0
        PB[PB<0] = 0

        # calculate x
        x = -V*PA*np.sqrt(PB)

        # calculate y
        y = np.ones_like(x)*float('nan')
        for i, t in enumerate(t_meas_block):
            # use the appropriate backward differences
            if t == 1:
                # this is the first data point in the experiment
                y[i] = (nA[i] - nA0[i])/t
            else:
                # this is not the first data point in the experiment
                y[i] = (nA[i] - nA[i-1])/(t - t_meas_block[i-1])

        # fit a straight line through the origin to the x-y data for this block
        param, param_ci, r_squared, y_pred = lls_parameters(y, x
                , model_has_intercept=False, use_rel_errors=False)
        
        # add the results to the model plot
        plt.plot(x, y, marker=markers[iBlock], ls='', color=colors[iBlock]
            , markerfacecolor='none', label=f'{blockT - 273.15:.0f} °C')
        plt.plot(x,y_pred, ls='-', color=colors[iBlock])
        
        # save the results for this block
        k[iBlock] = param[0]
        k_CI_lower[iBlock] = param_ci[0][0]
        k_CI_upper[iBlock] = param_ci[0][1]
        model_r_sq[iBlock] = r_squared
    
    # complete the model plots
    plt.xlabel('x (cm$^3$ atm$^{1.5}$)')
    plt.ylabel('y (mol min$^{-1}$)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('png/model_plot.png', dpi=300)
    plt.show(block=False)

    # tabulate, show, and save the model plot results
    results_df = pd.DataFrame({'T': block_T_values - 273.15, 'k': k
        , 'k_lower': k_CI_lower, 'k_upper': k_CI_upper, 'Rsq': model_r_sq})
    print('')
    print(results_df)
    print('')
    results_df.to_csv('model_plot_results.csv', index=False)

    # fit the Arrhenius expression to the T-k data
    k0, k0_CI, E, E_CI, r_sq_Arr = Arrhenius_parameters(k
            , block_T_values, Ren)

    # tabulate, show, and save the results
    data = [['k0', f'{k0:.3g}', f'{k0_CI[0]:.3g}', f'{k0_CI[1]:.3g}'
            , 'mol cm^-3^ min^-1^ atm^-1.5^'],
        ['E', f'{E:.3g}', f'{E_CI[0]:.3g}', f'{E_CI[1]:.3g}', 'kcal mol^-1^'],
        ['R_squared', f'{r_sq_Arr:.3g}', '','','']]
    result_df = pd.DataFrame(data, columns=['Parameter','value', 'lower_limit'
        , 'upper_limit','units'])
    print('')
    print(result_df)
    result_df.to_csv('arrhenius_results.csv', index=False)

    # generate, show, and save the Arrhenius plot
    k_pred = k0*np.exp(-E/Ren/block_T_values)
    plt.figure()
    plt.semilogy(1/block_T_values,k,color='tab:blue',marker='o'
        , markerfacecolor='none', ls='none')
    plt.semilogy(1/block_T_values,k_pred,color='k')
    plt.xlabel('T$^{-1}$ (K$^{-1}$)')
    plt.ylabel('k (min$^{-1}$)')
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.savefig('png/Arrhenius_plot.png', dpi=300)
    plt.show()

if __name__=="__main__":
    deliverables()
