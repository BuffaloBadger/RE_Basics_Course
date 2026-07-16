"""Post processing for Example 11.5.6 from REB, The Book"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from great_tables import GT, md, html

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# table with the first six rows of the experimental data
table_df = pd.read_csv('example_11_5_6_data.csv')
table_df = table_df.head(6)
data_tbl = (GT(table_df)
    .tab_stub(rowname_col='Experiment')
    .tab_stubhead(label='Experiment')
    .tab_spanner(
        label='Adjusted Inputs',
        columns = ['T', 'CA0', 't_meas']
    )
    .tab_spanner(
        label= 'Response',
        columns = ['CA_meas']
    )
    .fmt_number(columns='Experiment', decimals=0)
    .cols_label(T=html('T<br>(°C)'),
                CA0=html('C<sub>A,0</sub><br>(M)'),
                t_meas=html('t<sub>meas</sub><br>(min)'),
                CA_meas=html('C<sub>A,meas</sub><br>(M)')
    )
    .cols_align(
        align='center',
        columns=[0,1,2,3,4]
    )
)
data_tbl.gtsave('png/data_table.png', zoom=4.0)

# Arrhenius parameter table from fitting function
# read the data
df = pd.read_csv('results/ff_results.csv')
# drop the last column of units
df = df.drop('4', axis=1)
# rename the second and third columns
df = df.rename(columns = {'1': 'Value', '2': '95% CI'})
# add units to the first column and format as html
df['0'] = ['k<sub>0</sub> (min<sup>-1</sup>)', 'E (kJ mol<sup>-1</sup>)', 'R<sup>2</sup>']
# build the table
ff_Arr_tbl = (GT(df)
#    .tab_header(title='Arrhenius Parameters')
    .tab_stub(rowname_col='0')
    .tab_stubhead('Parameter')
    .fmt_scientific(rows=[0], exp_style='x10n')
    .fmt_number(rows=[1], decimals=1)
    .fmt_number(columns=[1], rows=[2], decimals=3)
    .cols_merge(columns=[2,3], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[2])
    .cols_align(align='center', columns=[1,2,3])
    .sub_missing(columns=[2,3], rows=[2], missing_text=' ')
)
ff_Arr_tbl.gtsave('png/ff_Arr_table.png', zoom=4.0)

# Arrhenius parameter table from linearized model
# read the data
df = pd.read_csv('results/lm_results.csv')
# drop the last column of units
df = df.drop('4', axis=1)
# rename the second and third columns
df = df.rename(columns = {'1': 'Value', '2': '95% CI'})
# add units to the first column and format as html
df['0'] = ['k<sub>0</sub> (min<sup>-1</sup>)'
        , 'E (kJ mol<sup>-1</sup>)'
        , 'R<sup>2</sup> (Arrhenius Plot)'
        , 'R<sup>2</sup> (Full Data Set)']
# build the table
lm_Arr_tbl = (GT(df)
#    .tab_header(title='Arrhenius Parameters')
    .tab_stub(rowname_col='0')
    .tab_stubhead('Parameter')
    .fmt_scientific(rows=[0], exp_style='x10n')
    .fmt_number(rows=[1], decimals=1)
    .fmt_number(columns=[1], rows=[2,3], decimals=3)
    .cols_merge(columns=[2,3], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[2,3])
    .cols_align(align='center', columns=[1,2,3])
    .sub_missing(columns=[2,3], rows=[2,3], missing_text=' ')
)
lm_Arr_tbl.gtsave('png/lm_Arr_table.png', zoom=4.0)

# Model plots parameter table from linearized model
# read the data
df = pd.read_csv('results/lm_model_plot_results.csv')
# build the table
lm_model_plot_tbl = (GT(df)
#    .tab_header(title='Model Plots Parameters')
    .fmt_number(columns=[0], decimals=0)
    .fmt_number(columns=[1,2,3,4], n_sigfig=3)
    .cols_label(block_T_values='T (°C)',
                lm_k=html('k (min<sup>-1</sup>)'),
                lm_k_CI_lower ='95% CI',
                lm_model_r_sq=html('R<sup>2</sup>'))
    .cols_merge(columns=[2,3], rows=[0,1,2,3]
                , pattern='[{0}, {1}]')
    .cols_align(align='center', columns=[0,1,2,3,4])
)
lm_model_plot_tbl.gtsave('png/lm_model_plots_table.png', zoom=4.0)

# Arrhenius parameter table from linearized approximate model
# read the data
df = pd.read_csv('results/am_results.csv')
# drop the last column of units
df = df.drop('4', axis=1)
# rename the second and third columns
df = df.rename(columns = {'1': 'Value', '2': '95% CI'})
# add units to the first column and format as html
df['0'] = ['k<sub>0</sub> (min<sup>-1</sup>)'
        , 'E (kJ mol<sup>-1</sup>)'
        , 'R<sup>2</sup> (Arrhenius Plot)'
        , 'R<sup>2</sup> (Full Data Set)']
# build the table
am_Arr_tbl = (GT(df)
#    .tab_header(title='Arrhenius Parameters')
    .tab_stub(rowname_col='0')
    .tab_stubhead('Parameter')
    .fmt_scientific(rows=[0], exp_style='x10n')
    .fmt_number(rows=[1], decimals=1)
    .fmt_number(columns=[1], rows=[2,3], decimals=3)
    .cols_merge(columns=[2,3], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[2,3])
    .cols_align(align='center', columns=[1,2,3])
    .sub_missing(columns=[2,3], rows=[2,3], missing_text=' ')
)
am_Arr_tbl.gtsave('png/am_Arr_table.png', zoom=4.0)

# Model plots parameter table from linearized approximate model
# read the data
df = pd.read_csv('results/am_model_plot_results.csv')
# build the table
am_model_plot_tbl = (GT(df)
#    .tab_header(title='Model Plots Parameters')
    .fmt_number(columns=[0], decimals=0)
    .fmt_number(columns=[1,2,3,4], n_sigfig=3)
    .cols_label(block_T_values='T (°C)',
                am_k=html('k (min<sup>-1</sup>)'),
                am_k_CI_lower ='95% CI',
                am_model_r_sq=html('R<sup>2</sup>'))
    .cols_merge(columns=[2,3], rows=[0,1,2,3]
                , pattern='[{0}, {1}]')
    .cols_align(align='center', columns=[0,1,2,3,4])
)
am_model_plot_tbl.gtsave('png/am_model_plots_table.png', zoom=4.0)

# combined model plot parameters table
# read the linear model data
df = pd.read_csv('results/lm_model_plot_results.csv')
# read the approximate model data
app_df = pd.read_csv('results/am_model_plot_results.csv')
# copy columns to the first dataframe
df['am_k'] = app_df['am_k']
df['am_k_CI_lower'] = app_df['am_k_CI_lower']
df['am_k_CI_upper'] = app_df['am_k_CI_upper']
df['am_model_r_sq'] = app_df['am_model_r_sq']
# build the table
comb_model_params_tbl = (GT(df)
#    .tab_header(title='Model Plots Parameter Comparison')
    .fmt_number(columns=[0], decimals=0)
    .fmt_number(columns=[1,2,3,4,5,6,7,8], n_sigfig=3)
    .cols_label(block_T_values='T (°C)',
                lm_k=html('k (min<sup>-1</sup>)'),
                lm_k_CI_lower ='95% CI',
                lm_model_r_sq=html('R<sup>2</sup>'),
                am_k=html('k (min<sup>-1</sup>)'),
                am_k_CI_lower ='95% CI',
                am_model_r_sq=html('R<sup>2</sup>'))
    .cols_merge(columns=[2,3], rows=[0,1,2,3]
                , pattern='[{0}, {1}]')
    .cols_merge(columns=[6,7], rows=[0,1,2,3]
                , pattern='[{0}, {1}]')
    .tab_spanner(label='Linear Model', columns=[1,2,3,4])
    .tab_spanner(label='Approximate Model', columns=[5,6,7,8])
    .cols_align(align='center', columns=[0,1,2,3,4,5,6,7,8])
)
comb_model_params_tbl.gtsave('png/combined_model_plot_params_table.png'
        , zoom=4.0)

# combined Arrhenius parameters
# read the linear model data
df = pd.read_csv('results/lm_results.csv')
# drop the last column of units
df = df.drop('4', axis=1)
# read the approximate model data
app_df = pd.read_csv('results/am_results.csv')
# copy columns to the first dataframe
df['4'] = app_df['1']
df['5'] = app_df['2']
df['6'] = app_df['3']
# add units to the first column and format as html
df['0'] = ['k<sub>0</sub> (min<sup>-1</sup>)'
        , 'E (kJ mol<sup>-1</sup>)'
        , 'R<sup>2</sup> (Arrhenius Plot)'
        , 'R<sup>2</sup> (Full Data Set)']
# rename the second, third, 4th abd 5tg=h columns
df = df.rename(columns = {'1': 'Value', '2': '95% CI', '4': 'v2', '5': 'v3'})
# build the table
comb_Arr_tbl = (GT(df)
#    .tab_header(title='Arrhenius Parameters')
    .tab_stub(rowname_col='0')
    .tab_stubhead('Parameter')
    .fmt_scientific(rows=[0], exp_style='x10n')
    .fmt_number(rows=[1], decimals=1)
    .fmt_number(columns=[1,4], rows=[2,3], decimals=3)
    .cols_merge(columns=[2,3], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[2,3])
    .cols_merge(columns=[5,6], rows=[0,1], pattern='[{0}, {1}]')
    .cols_merge(columns=[5,6], rows=[2,3])
    .cols_label(v2='Value', v3='95% CI')
    .cols_align(align='center', columns=[1,2,3,4,5,6])
    .tab_spanner(label='Linear Model', columns=[1,2,3])
    .tab_spanner(label='Approximate Model', columns=[4,5,6])
    .sub_missing(columns=[2,3,5,6], rows=[2,3], missing_text=' ')
)
comb_Arr_tbl.gtsave('png/combined_Arrhenius_plot_params.png', zoom=4.0)

# parameter estimation summary table
# read the data
ff_df = pd.read_csv('results/ff_results.csv')
lm_df = pd.read_csv('results/lm_results.csv')
am_df = pd.read_csv('results/am_results.csv')
# create the df for the table
df = pd.DataFrame()
df['Method'] = ['fitting function', 'linearized model', 'approximate model']
df['k0'] = [ff_df.iat[0,1], lm_df.iat[0,1], am_df.iat[0,1]]
df['k0_ll'] = [ff_df.iat[0,2], float('nan'), float('nan')]
df['k0_ul'] = [ff_df.iat[0,3], float('nan'), float('nan')]
df['E'] = [ff_df.iat[1,1], lm_df.iat[1,1], am_df.iat[1,1]]
df['E_ll'] = [ff_df.iat[1,2], float('nan'), float('nan')]
df['E_ul'] = [ff_df.iat[1,3], float('nan'), float('nan')]
df['R_sq'] = [ff_df.iat[2,1], lm_df.iat[3,1], am_df.iat[3,1]]
# build the table
summary_tbl = (GT(df)
#    .tab_header(title='Estimated Parameters Comparison')
    .fmt_scientific(columns=[1,2,3], exp_style='x10n')
    .fmt_number(columns=[4,5,6], decimals=1)
    .fmt_number(columns=[7], decimals=3)
    .cols_merge(columns=[2,3], rows=[0], pattern='[{0}, {1}]')
    .cols_merge(columns=[2,3], rows=[1,2])
    .cols_merge(columns=[5,6], rows=[0], pattern='[{0}, {1}]')
    .cols_merge(columns=[5,6], rows=[1,2])
    .sub_missing(columns=[2,3], rows=[1,2], missing_text=' ')
    .sub_missing(columns=[5,6], rows=[1,2], missing_text=' ')
    .cols_label(k0=html('k<sub>0</sub> (min<sup>-1</sup>)')
                ,k0_ll='95% CI', E=html('E (kJ mol<sup>-1</sup>)')
                ,E_ll='95% CI', R_sq=html('R<sup>2</sup>'))
    .cols_align(align='center', columns=[1,2,3,4,5,6,7])
)
summary_tbl.gtsave('png/param_est_summary.png', zoom=4.0)

# combined Arrhenius plots
lm_df = pd.read_csv('results/lm_Arr_plot_data.csv')
am_df = pd.read_csv('results/am_Arr_plot_data.csv')
T = lm_df['block_T_values'].to_numpy()
k_lm = lm_df['lm_k'].to_numpy()
k_pred_lm = lm_df['lm_k_pred'].to_numpy()
k_am = am_df['am_k'].to_numpy()
k_pred_am = am_df['am_k_pred'].to_numpy()
plt.figure()
plt.semilogy(1/T,k_lm,color='tab:blue',marker='o', ls='none'
        , label='linearize model')
plt.semilogy(1/T,k_pred_lm,color='tab:blue', ls='-')
plt.semilogy(1/T,k_am,color='tab:orange',marker='o', ls='none'
        , label='approximate model')
plt.semilogy(1/T,k_pred_am,color='tab:orange', ls='-')
plt.xlabel('T$^{-1}$ (K$^{-1}$)')
plt.ylabel('k (min$^{-1}$)')
plt.xticks(rotation=25)
plt.legend()
plt.tight_layout()
plt.savefig('png/combined_Arrhenius_plot.png', dpi=300)
plt.show(block=False)

# combined parity plots
df = pd.read_csv('example_11_5_6_data.csv')
CA_meas = df['CA_meas'].to_numpy()
T = df['T'].to_numpy()
CA0 = df['CA0'].to_numpy()
t_meas = df['t_meas'].to_numpy()
df = pd.read_csv('results/ff_plot_data.csv')
CA_ff = df['ff_CA_pred'].to_numpy()
eps_ff = df['ff_epsilon_expt']
df = pd.read_csv('results/lm_parity_plot_data.csv')
CA_lm = df['lm_CA_pred'].to_numpy()
eps_lm = df['lm_epsilon_expt'].to_numpy()
df = pd.read_csv('results/am_parity_plot_data.csv')
CA_am = df['am_CA_pred'].to_numpy()
eps_am = df['am_epsilon_expt'].to_numpy()
plt.figure()
plt.plot(CA_meas, CA_ff, color='tab:blue', marker='o', ls=''
        , label='fitting function', markerfacecolor='none')
plt.plot(CA_meas, CA_lm, color='tab:orange', marker='x', ls=''
        , label='linearized model')
plt.plot(CA_meas, CA_am, color='tab:green', marker='+', ls=''
        , label='approximate model')
plt.plot([min(CA_meas),max(CA_meas)],[min(CA_meas),max(CA_meas)]
        , color = 'k', ls = '-', label = 'Parity Line')
plt.xlabel('$C_{A,meas}$ (M)')
plt.ylabel('$C_{A,pred}$ (M)')
plt.legend()
plt.savefig('./png/combined_parity_plot.png', dpi=300)
plt.show(block=False)

# combined residuals plots
plt.figure()
plt.plot(CA0, eps_ff, markerfacecolor='none', marker='o'
            , color='tab:blue', ls='', label='fitting function')
plt.plot(CA0, eps_lm, markerfacecolor='none', marker='x'
            , color='tab:orange', ls='', label='linearized model')
plt.plot(CA0, eps_am, markerfacecolor='none', marker='+'
            , color='tab:green', ls='', label='approximate model')
plt.axhline(y=0, color='k')
plt.xlabel('$C_{A,0}$ (M)')
plt.ylabel('Residual (M)')
plt.legend()
plt.tight_layout()
plt.savefig('./png/combined_CA_residuals_plot.png', dpi=300)
plt.show(block=False)

plt.figure()
plt.plot(T, eps_ff, markerfacecolor='none', marker='o'
            , color='tab:blue', ls='', label='fitting function')
plt.plot(T, eps_lm, markerfacecolor='none', marker='x'
            , color='tab:orange', ls='', label='linearized model')
plt.plot(T, eps_am, markerfacecolor='none', marker='+'
            , color='tab:green', ls='', label='approximate model')
plt.axhline(y=0, color='k')
plt.xlabel('T (°C)')
plt.ylabel('Residual (M)')
plt.legend()
plt.tight_layout()
plt.savefig('./png/combined_T_residuals_plot.png', dpi=300)
plt.show(block=False)

plt.figure()
plt.plot(t_meas, eps_ff, markerfacecolor='none', marker='o'
            , color='tab:blue', ls='', label='fitting function')
plt.plot(t_meas, eps_lm, markerfacecolor='none', marker='x'
            , color='tab:orange', ls='', label='linearized model')
plt.plot(t_meas, eps_am, markerfacecolor='none', marker='+'
            , color='tab:green', ls='', label='approximate model')
plt.axhline(y=0, color='k')
plt.xlabel('$t_{meas}$ (min)')
plt.ylabel('Residual (M)')
plt.legend()
plt.tight_layout()
plt.savefig('./png/combined_time_residuals_plot.png', dpi=300)
plt.show()

'''

    # combined parity plots
    plt.figure()
    plt.plot(CA_meas, ff_CA_pred, marker='x', ls='', color='tab:blue'
             ,markerfacecolor='none', label='fitting function')
    plt.plot(CA_meas, lm_CA_pred, marker='+', ls='', color='tab:orange'
             , markerfacecolor='none', label='linear model')
    plt.plot(CA_meas, am_CA_pred, marker='o', ls='', color='tab:green'
             , markerfacecolor='none', label='approximate model')
    plt.plot([min(CA_meas),max(CA_meas)],[min(CA_meas),max(CA_meas)]
             , color = 'k', ls = '-', label = 'Parity Line')
    plt.xlabel('$C_{A,meas}$ (M)')
    plt.ylabel('$C_{A,pred}$ (M)')
    plt.legend()
    plt.savefig('./pdf/example_11_5_6_combined_parity.pdf')
    plt.savefig('./png/example_11_5_6_combined_parity.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_6_combined_parity.png', dpi=300)
    plt.show(block=False)

    # combined residuals plots
    plt.figure()
    plt.plot(CA_0, ff_epsilon_expt, markerfacecolor='none', marker='x'
             , color='tab:blue', ls='', label='fitting function')
    plt.plot(CA_0, lm_epsilon_expt, markerfacecolor='none', marker='+'
             , color='tab:orange', ls='', label='linear model')
    plt.plot(CA_0, epsilon_expt_am, markerfacecolor='none', marker='o'
             , color='tab:green', ls='', label='approximate model')
    plt.axhline(y=0, color='k')
    plt.xlabel('$C_{A,0}$ (M)')
    plt.ylabel('Residual (M)')
    plt.legend()
    plt.savefig('./pdf/example_11_5_6_CA_residuals.pdf')
    plt.savefig('./png/example_11_5_6_CA_residuals.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_6_CA_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure()
    plt.plot(T-273.15, ff_epsilon_expt, markerfacecolor='none', marker='x'
             , color='tab:blue', ls='', label='fitting function')
    plt.plot(T-273.15, lm_epsilon_expt, markerfacecolor='none', marker='+'
             , color='tab:orange', ls='', label='linear model')
    plt.plot(T-273.15, epsilon_expt_am, markerfacecolor='none', marker='o'
             , color='tab:green', ls='', label='approximate model')
    plt.axhline(y=0, color='k')
    plt.xlabel('T (°C)')
    plt.ylabel('Residual (M)')
    plt.legend()
    plt.savefig('./pdf/example_11_5_6_T_residuals.pdf')
    plt.savefig('./png/example_11_5_6_T_residuals.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_6_T_residuals.png', dpi=300)
    plt.show(block=False)

    plt.figure()
    plt.plot(t_meas, ff_epsilon_expt, markerfacecolor='none', marker='x'
             , color='tab:blue', ls='', label='fitting function')
    plt.plot(t_meas, lm_epsilon_expt, markerfacecolor='none', marker='+'
             , color='tab:orange', ls='', label='linear model')
    plt.plot(t_meas, epsilon_expt_am, markerfacecolor='none', marker='o'
             , color='tab:green', ls='', label='approximate model')
    plt.axhline(y=0, color='k')
    plt.xlabel('$t_{meas}$ (min)')
    plt.ylabel('Residual (M)')
    plt.legend()
    plt.savefig('./pdf/example_11_5_6_time_residuals.pdf')
    plt.savefig('./png/example_11_5_6_time_residuals.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_6_time_residuals.png', dpi=300)
    plt.show(block=False)
    
    # combined model plots for linear model
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
    plt.savefig('./pdf/example_11_5_6_linear_model_plots.pdf')
    plt.savefig('./png/example_11_5_6_linear_model_plots.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_6_linear_model_plots.png', dpi=300)
    plt.show(block=False)

    # combined model plots for approximate model
    plot_data = np.transpose(np.array(am_model_plot_data))
    plt.figure()
    colors = ['tab:blue','tab:orange','tab:green','tab:purple']
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
    plt.legend(title='approximate model')
    plt.savefig('./pdf/example_11_5_6_approx_model_plots.pdf')
    plt.savefig('./png/example_11_5_6_approx_model_plots.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_6_approx_model_plots.png', dpi=300)
    plt.show(block=False)



    # combined Arrhenius plots
    T_K = block_T_values + 273.15
    lm_k_pred = lm_k0*np.exp(-lm_E/(R*T_K))
    am_k_pred = am_k0*np.exp(-am_E/(R*T_K))
    plt.figure()
    plt.semilogy(1/T_K, lm_k, color='tab:blue', marker='o', ls='none'
            , label=f'linear model, R$^2$ = {lm_r_sq_Arr:.3f}')
    plt.semilogy(1/T_K, lm_k_pred, color='tab:blue', ls='-')
    plt.semilogy(1/T_K, am_k, color='tab:orange', marker='o', ls='none'
            , label=f'approximate model, R$^2$ = {am_r_sq_Arr:.3f}')
    plt.semilogy(1/T_K, am_k_pred, color='tab:orange', ls='-')
    plt.xlabel('T$^{-1}$ (K$^{-1}$)')
    plt.ylabel('k (min$^{-1}$)')
    plt.xticks(rotation=25)
    plt.legend()
    plt.tight_layout()
    plt.savefig('./pdf/example_11_5_6_Arrhenius_plots.pdf')
    plt.savefig('./png/example_11_5_6_Arrhenius_plots.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_6_Arrhenius_plots.png', dpi=300)
    plt.show()
    
'''