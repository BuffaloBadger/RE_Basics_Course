"""Calculations for the Class 27 Learning Activity in REB, The Course"""

# import libraries
import pandas as pd
import numpy as np
from reb_utils import lls_parameters
from reb_utils import Arrhenius_parameters
import matplotlib.pyplot as plt
import os

# global constants available to all functions
# given
V = 1.0 # gal
# known
R = 1.987 # BTU lbmol^-1^ °R^-1^

# read the data file
expt_df = pd.read_csv('activity_27_data.csv')

# split it into three dataframes, one for each measurement time
df_10 = expt_df[['T', 'CA0', 'CB0', 'CZ_10']].copy()
df_10.rename(columns={'CZ_10' : 'CZ_meas'}, inplace=True)
df_30 = expt_df[['T', 'CA0', 'CB0', 'CZ_30']].copy()
df_30.rename(columns={'CZ_30' : 'CZ_meas'}, inplace=True)
df_50 = expt_df[['T', 'CA0', 'CB0', 'CZ_50']].copy()
df_50.rename(columns={'CZ_50' : 'CZ_meas'}, inplace=True)

# add columns containing the measurement time
t10 = np.ones(len(df_10)) * 10.0
df_10['t_meas'] = t10.tolist()
t30 = np.ones(len(df_30)) * 30.0
df_30['t_meas'] = t30.tolist()
t50 = np.ones(len(df_50)) * 50.0
df_50['t_meas'] = t50.tolist()

# recombine into a single dataframe
expt_df = pd.concat([df_10, df_30, df_50])

# extract the data as vectors
T = expt_df['T'].to_numpy() # °C
T = (T + 273.15)*1.8 # °R
CA0 = expt_df['CA0'].to_numpy()
CB0 = expt_df['CB0'].to_numpy()
t_meas = expt_df['t_meas'].to_numpy()
CZ_meas = expt_df['CZ_meas'].to_numpy()

# deliverables function
def deliverables():
    # make sure a png folder exists
    if not os.path.isdir('png'):
        os.makedirs('./png')

    # get the block temperatures
    block_T_values = np.array(expt_df['T'].unique())

    # allocate storage for the model plots results
    k = np.ones_like(block_T_values)*float('nan')
    k_CI_lower = np.ones_like(block_T_values)*float('nan')
    k_CI_upper = np.ones_like(block_T_values)*float('nan')
    model_r_sq = np.ones_like(block_T_values)*float('nan')
    model_plot_data = []

    # loop through the data blocks
    for iBlock, T_C in enumerate(block_T_values):
        # extract the same-temperature data block
        block_df = expt_df[expt_df['T'] == T_C]

        # extract the data in the block
        CA0_block = block_df['CA0'].to_numpy()
        CB0_block = block_df['CB0'].to_numpy()
        t_meas_block = block_df['t_meas'].to_numpy()
        CZ_meas_block = block_df['CZ_meas'].to_numpy()

        # calculate x and add it to the data for this block
        x = -t_meas_block
        model_plot_data.append(x)

        # calculate y and add it to the data for this block
        y = np.ones_like(CZ_meas_block)*float('nan')
        for i, CZ in enumerate(CZ_meas_block):
            CA = CA0[i] - + CZ
            if CA0_block[i] == CB0_block[i]:
                y[i] = 1/CA0_block[i] - 1/CA
            else:
                y[i] = 1/(CA0_block[i] - CB0_block[i]) \
                    * np.log(CA0_block[i] *(CB0_block[i] - CA0_block[i] + CA)
                    / (CB0_block[i] * CA))
        model_plot_data.append(y)

        # fit a straight line through the origin to the x-y data for this block
        param, param_ci, r_squared, y_pred = lls_parameters(y, x
                , model_has_intercept=False, use_rel_errors=False)
        
        # add y_pred to the data for this block
        model_plot_data.append(y_pred)
        
        # save the results for this block
        k[iBlock] = param[0]
        k_CI_lower[iBlock] = param_ci[0][0]
        k_CI_upper[iBlock] = param_ci[0][1]
        model_r_sq[iBlock] = r_squared

    # generate, show and save a model plots results table
    df = pd.DataFrame({'T' : block_T_values
                       ,'k' : k
                       ,'k_lower' : k_CI_lower
                       ,'k_upper' : k_CI_upper
                       ,'R_sq' : model_r_sq})
    print('')
    print(df)
    df.to_csv('model_plot_results.csv', index=False)

    # generate, show, and save a combined model plot
    plot_data = np.transpose(np.array(model_plot_data))
    plt_colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple']
    plt_markers = ['x', 'o', '+', '^']
    plt.figure()
    count = 0
    for i, T_C in enumerate(block_T_values):
        plt.plot(plot_data[:,3*count], plot_data[:,3*count+1]
            , marker=plt_markers[count], ls='', color=plt_colors[count]
            , label=f'{T_C:.0f} °C', markerfacecolor='none')
        plt.plot(plot_data[:,3*count], plot_data[:,3*count+2]
        , ls='-', color=plt_colors[count])
        count = count + 1
    plt.xlabel('x (min)')
    plt.ylabel('y (gal lbmol$^{-1}$)')
    plt.legend()
    plt.savefig('png/model_plots.png', dpi=300)
    plt.show(block=False)

    # fit the Arrhenius expression to the T-K data
    k0, k0_CI, E, E_CI, r_sq_Arr = Arrhenius_parameters(k
        , (block_T_values + 273.15)*1.8, R)

    # generate, show and save an Arrhenius results table
    results = [
        ['k0', k0, k0_CI[0], k0_CI[1], 'gal/lbmol/min']
        , ['E', E, E_CI[0], E_CI[1], 'BTU/lbmol']
        , ['R^2 Arr', r_sq_Arr, float('nan'), float('nan'),'']
    ]
    df = pd.DataFrame(results, columns=['Parameter', 'Value', 'lower_limit'
            , 'upper_limit', 'units'])
    df.to_csv('Arrhenius_results.csv', index=False)
    print('')
    print(df)

    # generate, show, and save an Arrhenius plot
    T_R = (block_T_values + 273.15)*1.8
    k_pred = k0*np.exp(-E/(R*T_R))
    plt.figure()
    plt.semilogy(1/block_T_values,k,color='tab:blue',marker='o', ls='none')
    plt.semilogy(1/block_T_values,k_pred,color='k')
    plt.xlabel('T$^{-1}$ (°R$^{-1}$)')
    plt.ylabel('k (gal lbmol$^{-1}$ min$^{-1}$)')
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.savefig('png/Arrhenius_plot.png', dpi=300)
    plt.show(block=False)

    # calculate the predicted responses
    k = k0*np.exp(-E/(R*T))
    x = -t_meas
    CZ_pred = np.ones_like(CZ_meas)*float('nan')
    for i, CZ in enumerate(CZ_meas):
        if CA0[i] == CB0[i]:
            CZ_pred[i] = CA0[i] - 1/(1/CA0[i] - k[i]*x[i])
        else:
            CZ_pred[i] = CA0[i] - (CA0[i] - CB0[i])/(1 - CB0[i]/CA0[i]
                * np.exp((CA0[i] - CB0[i])*k[i]*x[i]))
    
    # calculate the  coefficient of determination
    CZ_mean = np.mean(CZ_meas)
    ss_res = np.sum(np.square(CZ_meas - CZ_pred))
    ss_tot = np.sum(np.square(CZ_meas - CZ_mean))
    r_sq = 1 - ss_res/ss_tot

    # generate, show and save a parity plot
    plt.figure()
    plt.plot(CZ_meas, CZ_pred, marker='x', ls='', color='tab:blue'
             ,markerfacecolor='none', label='Data')
    plt.plot([min(CZ_meas),max(CZ_meas)],[min(CZ_meas),max(CZ_meas)]
             , color = 'k', ls = '-', label = 'Parity Line')
    plt.xlabel('$C_{Z,meas}$ (lbmol gal$^{-1}$)')
    plt.ylabel('$C_{Z,pred}$ (lbmol gal$^{-1}$)')
    plt.legend(title=f'R$^2$ = {r_sq:.3f}')
    plt.savefig('./png/parity.png', dpi=300)
    plt.show(block=False)

    # calculate the experiment residuals
    epsilon_expt = CZ_meas - CZ_pred

    # generate, show, and save residuals plots
    plt.figure() 
    plt.plot(CA0, epsilon_expt, color = 'tab:blue', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Initial Concentration of A (lbmol gal$^{-1}$)")
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.tight_layout()
    plt.savefig('png/CA_residuals.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(CB0, epsilon_expt, color = 'tab:blue', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Initial Concentration of B (lbmol gal$^{-1}$)")
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.tight_layout()
    plt.savefig('png/CB_residuals.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(T/1.8 - 273.15, epsilon_expt, color = 'tab:blue', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.tight_layout()
    plt.savefig('png/T_residuals.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(t_meas, epsilon_expt, color = 'tab:blue', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Measurement Time (min)")
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.tight_layout()
    plt.savefig('png/time_residuals.png')
    plt.show()

if __name__=="__main__":
    deliverables()
