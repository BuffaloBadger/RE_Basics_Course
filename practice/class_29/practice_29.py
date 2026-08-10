"""Calculations for the Class 29 Practice Assignment from REB, The Course"""

# import libraries
import pandas as pd
import numpy as np
from reb_utils import lls_parameters
from reb_utils import Arrhenius_parameters
import matplotlib.pyplot as plt
import os

# global constants available to all functions
# given
V = 1.0 # L (basis)
# known
Ren = 1.987E-3 # kcal/mol

# experimental data as arrays
expt_df = pd.read_csv('practice_29_data.csv')
t_meas_1 = expt_df['t_meas_exp1'].to_numpy()
CY_meas_1 = expt_df['CY_meas_exp1'].to_numpy()
t_meas_2 = expt_df['t_meas_exp2'].to_numpy()
CY_meas_2 = expt_df['CY_meas_exp2'].to_numpy()
t_meas_3 = expt_df['t_meas_exp3'].to_numpy()
CY_meas_3 = expt_df['CY_meas_exp3'].to_numpy()

# deliverables function
def deliverables():
    # make sure a png folder exists
    if not os.path.isdir('png'):
        os.makedirs('./png')

    # analyze the first experiment
    CA0 = 1.0
    CB0 = 0.9
    CA = CA0 - CY_meas_1
    CB = CB0 - CY_meas_1
    x = CA*CB
    y = np.ones_like(x)
    for i, CY in enumerate(CY_meas_1):
        if i == 0:
            y[i] = CY/t_meas_1[0]
        else:
            y[i] = (CY - CY_meas_1[i-1])/(t_meas_1[i] - t_meas_1[i-1])

    # fit a straight line through the origin to the x-y data for this block
    param, param_ci, r_squared, y_pred = lls_parameters(y, x
            , model_has_intercept=False, use_rel_errors=False)
    
    # extract the results
    k1 = param[0]
    k1_lower = param_ci[0][0]
    k1_upper = param_ci[0][1]
    rSq1 = r_squared

    # generate, show, and save the model plot
    plt.figure()
    plt.plot(x, y, marker='x', ls='', color='k')
    plt.plot(x,y_pred, ls='-', color='k')
    plt.xlabel('x (mol$^2$ L$^{-2}$)')
    plt.ylabel('y (mol L$^{-1}$ min$^{-1}$)')
    plt.tight_layout()
    plt.savefig('png/expt1_model_plot.png', dpi=300)
    plt.show(block=False)

    # analyze the second experiment
    CA0 = 1.0
    CB0 = 0.75
    CA = CA0 - CY_meas_2
    CB = CB0 - CY_meas_2
    x = CA*CB
    y = np.ones_like(x)
    for i, CY in enumerate(CY_meas_2):
        if i == 0:
            y[i] = CY/t_meas_2[0]
        else:
            y[i] = (CY - CY_meas_2[i-1])/(t_meas_2[i] - t_meas_2[i-1])

    # fit a straight line through the origin to the x-y data for this block
    param, param_ci, r_squared, y_pred = lls_parameters(y, x
            , model_has_intercept=False, use_rel_errors=False)
    
    # extract the results
    k2 = param[0]
    k2_lower = param_ci[0][0]
    k2_upper = param_ci[0][1]
    rSq2 = r_squared

    # generate, show, and save the model plot
    plt.figure()
    plt.plot(x, y, marker='x', ls='', color='k')
    plt.plot(x,y_pred, ls='-', color='k')
    plt.xlabel('x (mol$^2$ L$^{-2}$)')
    plt.ylabel('y (mol L$^{-1}$ min$^{-1}$)')
    plt.tight_layout()
    plt.savefig('png/expt2_model_plot.png', dpi=300)
    plt.show(block=False)

    # analyze the third experiment
    CA0 = 0.75
    CB0 = 1.0
    CA = CA0 - CY_meas_3
    CB = CB0 - CY_meas_3
    x = CA*CB
    y = np.ones_like(x)
    for i, CY in enumerate(CY_meas_3):
        if i == 0:
            y[i] = CY/t_meas_3[0]
        else:
            y[i] = (CY - CY_meas_3[i-1])/(t_meas_3[i] - t_meas_3[i-1])

    # fit a straight line through the origin to the x-y data for this block
    param, param_ci, r_squared, y_pred = lls_parameters(y, x
            , model_has_intercept=False, use_rel_errors=False)
    
    # extract the results
    k3 = param[0]
    k3_lower = param_ci[0][0]
    k3_upper = param_ci[0][1]
    rSq3 = r_squared

    # generate, show, and save the model plot
    plt.figure()
    plt.plot(x, y, marker='x', ls='', color='k')
    plt.plot(x,y_pred, ls='-', color='k')
    plt.xlabel('x (mol$^2$ L$^{-2}$)')
    plt.ylabel('y (mol L$^{-1}$ min$^{-1}$)')
    plt.tight_layout()
    plt.savefig('png/expt3_model_plot.png', dpi=300)
    plt.show(block=False)

    # tabulate, show, and save the model plot results
    results = [[70.0, k1, k1_lower, k1_upper, rSq1],
               [80.0, k2, k2_lower, k2_upper, rSq2],
               [90.0, k3, k3_lower, k3_upper, rSq3]]
    results_df = pd.DataFrame(results, columns = ['T', 'k', 'k_lower', 'k_upper'
            , 'Rsq'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('model_plot_results.csv', index=False)

    # fit the Arrhenius expression to the T-k data
    T = np.array([70.0, 80.0, 90.0]) + 273.15
    k = np.array([k1, k2, k3])
    k0, k0_CI, E, E_CI, r_sq_Arr = Arrhenius_parameters(k
            , T, Ren)

    # tabulate, show, and save the results
    data = [['k0', f'{k0:.3g}', f'{k0_CI[0]:.3g}', f'{k0_CI[1]:.3g}'
            , 'L mol^-1^ min^-1^'],
        ['E', f'{E:.3g}', f'{E_CI[0]:.3g}', f'{E_CI[1]:.3g}', 'kcal mol^-1^'],
        ['R_squared', f'{r_sq_Arr:.3g}', '','','']]
    result_df = pd.DataFrame(data, columns=['Parameter','value', 'lower_limit'
        , 'upper_limit','units'])
    print('')
    print(result_df)
    result_df.to_csv('arrhenius_results.csv', index=False)

    # generate, show, and save the Arrhenius plot
    k_pred = k0*np.exp(-E/Ren/T)
    plt.figure()
    plt.semilogy(1/T,k,color='k',marker='o', ls='none')
    plt.semilogy(1/T,k_pred,color='k')
    plt.xlabel('T$^{-1}$ (K$^{-1}$)')
    plt.ylabel('k (L mol$^{-1}$ min$^{-1}$)')
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.savefig('png/Arrhenius_plot.png', dpi=300)
    plt.show()

if __name__=="__main__":
    deliverables()
