"""Calculations for the Class 27 Practice Assignment in REB, The Course"""

# import libraries
import pandas as pd
import numpy as np
from reb_utils import lls_parameters
from reb_utils import Arrhenius_parameters
import matplotlib.pyplot as plt
import os

# global constants available to all functions
# given
V = 1.0 # L
# known
R = 1.987E-3 # kcal mol^-1^ K^-1^

# experimental data as arrays
expt_df = pd.read_csv('practice_27_data.csv')
T = expt_df['T'].to_numpy() + 273.15
CA0 = expt_df['CA0'].to_numpy()
CB0 = expt_df['CB0'].to_numpy()
t_meas = expt_df['t_meas'].to_numpy()
CA_meas = expt_df['CA_meas'].to_numpy()

# deliverables function
def deliverables():
    # make sure a png folder exists
    if not os.path.isdir('png'):
        os.makedirs('./png')

    # get the block temperatures
    block_T_values = np.array(expt_df['T'].unique()) + 273.15

    # allocate storage for the model plots results
    k = np.ones_like(block_T_values)*float('nan')
    k_CI_lower = np.ones_like(block_T_values)*float('nan')
    k_CI_upper = np.ones_like(block_T_values)*float('nan')
    model_r_sq = np.ones_like(block_T_values)*float('nan')
    model_plot_data = []

    # loop through the data blocks
    for iBlock, blockT in enumerate(block_T_values):
        # extract the same-temperature data block
        block_df = expt_df[expt_df['T'] + 273.15 == blockT]

        # extract the data in the block
        CA0_block = block_df['CA0'].to_numpy()
        CB0_block = block_df['CB0'].to_numpy()
        t_meas_block = block_df['t_meas'].to_numpy()
        CA_meas_block = block_df['CA_meas'].to_numpy()

        # calculate x and add it to the data for this block
        x = -t_meas_block
        model_plot_data.append(x)

        # calculate y and add it to the data for this block
        y = np.ones_like(CA_meas_block)*float('nan')
        for i, CA in enumerate(CA_meas_block):
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
    df = pd.DataFrame({'T' : block_T_values - 273.15
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
    for i, T_K in enumerate(block_T_values):
        plt.plot(plot_data[:,3*count], plot_data[:,3*count+1]
            , marker=plt_markers[count], ls='', color=plt_colors[count]
            , label=f'{T_K-273.15:.0f} °C', markerfacecolor='none')
        plt.plot(plot_data[:,3*count], plot_data[:,3*count+2]
        , ls='-', color=plt_colors[count])
        count = count + 1
    plt.xlabel('x (min)')
    plt.ylabel('y (L mol$^{-1}$)')
    plt.legend()
    plt.savefig('png/model_plots.png', dpi=300)
    plt.show(block=False)

    # fit the Arrhenius expression to the T-K data
    k0, k0_CI, E, E_CI, r_sq_Arr = Arrhenius_parameters(k
        , block_T_values, R)

    # generate, show and save an Arrhenius results table
    results = [
        ['k0', k0, k0_CI[0], k0_CI[1], 'L/mol/min']
        , ['E', E, E_CI[0], E_CI[1], 'kcal/mol']
        , ['R^2 Arr', r_sq_Arr, float('nan'), float('nan'),'']
    ]
    df = pd.DataFrame(results, columns=['Parameter', 'Value', 'lower_limit'
            , 'upper_limit', 'units'])
    df.to_csv('Arrhenius_results.csv', index=False)
    print('')
    print(df)

    # generate, show, and save an Arrhenius plot
    k_pred = k0*np.exp(-E/R/block_T_values)
    plt.figure()
    plt.semilogy(1/block_T_values,k,color='tab:blue',marker='o', ls='none')
    plt.semilogy(1/block_T_values,k_pred,color='k')
    plt.xlabel('T$^{-1}$ (K$^{-1}$)')
    plt.ylabel('k (L mol$^{-1}$ min$^{-1}$)')
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.savefig('png/Arrhenius_plot.png', dpi=300)
    plt.show(block=False)

    # calculate the predicted responses
    k = k0*np.exp(-E/(R*T))
    x = -t_meas
    CA_pred = np.ones_like(CA_meas)*float('nan')
    for i, CA in enumerate(CA_meas):
        if CA0[i] == CB0[i]:
            CA_pred[i] = 1/(1/CA0[i] - k[i]*x[i])
        else:
            CA_pred[i] = (CA0[i] - CB0[i])/(1 - CB0[i]/CA0[i]
                * np.exp((CA0[i] - CB0[i])*k[i]*x[i]))
    
    # calculate the  coefficient of determination
    CA_mean = np.mean(CA_meas)
    ss_res = np.sum(np.square(CA_meas - CA_pred))
    ss_tot = np.sum(np.square(CA_meas - CA_mean))
    r_sq = 1 - ss_res/ss_tot

    # generate, show and save a parity plot
    plt.figure()
    plt.plot(CA_meas, CA_pred, marker='x', ls='', color='tab:blue'
             ,markerfacecolor='none', label='Data')
    plt.plot([min(CA_meas),max(CA_meas)],[min(CA_meas),max(CA_meas)]
             , color = 'k', ls = '-', label = 'Parity Line')
    plt.xlabel('$C_{A,meas}$ (M)')
    plt.ylabel('$C_{A,pred}$ (M)')
    plt.legend(title=f'R$^2$ = {r_sq:.3f}')
    plt.savefig('./png/parity.png', dpi=300)
    plt.show(block=False)

    # calculate the experiment residuals
    epsilon_expt = CA_meas - CA_pred

    # generate, show, and save residuals plots
    plt.figure() 
    plt.plot(CA0, epsilon_expt, color = 'tab:blue', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Initial Concentration of A (mol L$^{-1}$)")
    plt.ylabel("Residual (mol L$^{-1}$)")
    plt.savefig('png/CA_residuals.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(CB0, epsilon_expt, color = 'tab:blue', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Initial Concentration of B (mol L$^{-1}$)")
    plt.ylabel("Residual (mol L$^{-1}$)")
    plt.savefig('png/CB_residuals.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(T-273.15, epsilon_expt, color = 'tab:blue', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Residual (mol L$^{-1}$)")
    plt.savefig('png/T_residuals.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(t_meas, epsilon_expt, color = 'tab:blue', marker='o'
            , markerfacecolor='none', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("Measurement Time (min)")
    plt.ylabel("Residual (mol L$^{-1}$)")
    plt.savefig('png/time_residuals.png')
    plt.show()

if __name__=="__main__":
    deliverables()
