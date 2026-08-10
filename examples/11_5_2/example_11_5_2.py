"""Calculations for Example 11.5.2 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ates
from reb_utils import fit_to_SR_data
import os

# global constants available to all functions
# given
V = 0.1 # L
# known
R = 1.987E-3 # kcal/mol/K

# experimental data as arrays
expt_df = pd.read_csv('example_11_5_2_data.csv')
T_K = expt_df['T'].to_numpy() # K
Vdot = (expt_df['Vdot'].to_numpy())*1.0E-3 # L /min
CA_in = expt_df['CA_in'].to_numpy() # mol /L
CB_in = expt_df['CB_in'].to_numpy() # mol /L
CY_in = expt_df['CY_in'].to_numpy() # mol /L
CZ_in = expt_df['CZ_in'].to_numpy() # mol /L
CY_out = expt_df['CY_out'].to_numpy() # mol /L

# global variables for the current values of k0 and E and the index of the 
# current experiment
g_k0 = float('nan')
g_E = float('nan')
g_iExpt = -1

# CSTR model function
def cstr_model_variables(iExpt):
    # make the current experiment's index available to the residuals function
    global g_iExpt
    g_iExpt = iExpt

    # define an initial guess for the CSTR model variables
    nA_in = CA_in[iExpt] * Vdot[iExpt]
    nB_in = CB_in[iExpt] * Vdot[iExpt]
    nY_in = CY_in[iExpt] * Vdot[iExpt]
    nZ_in = CZ_in[iExpt] * Vdot[iExpt]
    initial_guess = np.array([nA_in, nB_in, nY_in, nZ_in])

    # solve the cstr design equations
    soln, success, message = solve_ates(cstr_residuals, initial_guess)

    # check for solver issues
    if not success:
        print('')
        print(f'CSTR Model issue for experiment {iExpt}: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit')
    
    # return the CSTR model variables
    return soln[0], soln[1], soln[2], soln[3]

# CSTR residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA_out = guess[0]
    nB_out = guess[1]
    nY_out = guess[2]
    nZ_out = guess[3]

    # calculate the additional unknowns
    nA_in = CA_in[g_iExpt] * Vdot[g_iExpt]
    nB_in = CB_in[g_iExpt] * Vdot[g_iExpt]
    nY_in = CY_in[g_iExpt] * Vdot[g_iExpt]
    nZ_in = CZ_in[g_iExpt] * Vdot[g_iExpt]  
    k = g_k0 * np.exp(-g_E/(R*T_K[g_iExpt]))  
    CA = nA_out/Vdot[g_iExpt]
    CB = nB_out/Vdot[g_iExpt]
    r = k*CA*CB

    # evaluate the CSTR design equation residuals
    epsilon_1 = nA_in - nA_out - r*V
    epsilon_2 = nB_in - nB_out - r*V
    epsilon_3 = nY_in - nY_out + r*V
    epsilon_4 = nZ_in - nZ_out + r*V

    return np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4])

# predicted responses function
def predicted_responses(adj_inputs, param, E):
    # make k0 and E available to the residuals function
    global g_k0, g_E
    g_k0 = 10**param
    g_E = E

    # allocate storage for the predicted responses
    CY_out_pred = np.ones_like(T_K)*float('nan')

    # loop through the experiments
    for iExpt, T in enumerate(T_K):
        # solve the CSTR design equations
        nA_out, nB_out, nY_out, nZ_out = cstr_model_variables(iExpt)

        # calculate the predicted response
        CY_out_pred[iExpt] = nY_out/Vdot[g_iExpt]
    
    # return the predicted responses
    return CY_out_pred

# deliverables function
def deliverables():
    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([T_K, Vdot, CA_in, CB_in, CY_in, CZ_in]))

    # define a guess for the base 10 log of k0 and E
    par_guess = [0.0, 10.0]

    # estimate the parameters
    param, param_ci, r_squared, CY_out_pred = fit_to_SR_data(par_guess
        , adj_inputs, CY_out, predicted_responses, use_rel_error=False)
    
    # extract the parameter estimates and their confidence intervals
    k0 = 10**param[0]
    k0_CI = 10**param_ci[0,:]
    E = param[1]
    E_CI = param_ci[1,:]

    # calculate the experiment residuals
    epsilon_expt = CY_out - CY_out_pred

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
    reb_book_path = '../../../RE_Basics/solutions/ch11_ex2/'
    if not os.path.isdir(reb_book_path):
        # create the folder
        os.makedirs(reb_book_path)

    # tabulate, show and save the fitting results
    fitting_results = [['k0', f'{k0:.3g}', 'L mol^-1^ min^-1^'],
        ['k_lower_limit', f'{k0_CI[0]:.3g}', 'L mol^-1^ min^-1^'],
        ['k_upper_limit', f'{k0_CI[1]:.3g}', 'L mol^-1^ min^-1^'],
        ['E', f'{E:.3g}', 'kcal mol^-1^'],
        ['E_lower_limit', f'{E_CI[0]:.3g}', 'kcal mol^-1^'],
        ['E_upper_limit', f'{E_CI[1]:.3g}', 'kcal mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    fitting_results_df = pd.DataFrame(fitting_results, columns=['quantity'
            ,'value','units'])
    print(" ")
    print(fitting_results_df)
    print(" ")
    fitting_results_df.to_csv("./csv/example_11_5_2_results.csv", index=False)
    fitting_results_df.to_csv(reb_book_path + 
        "example_11_5_2_results.csv", index=False)

    # generate, show and save a parity plot
    plt.figure(1) 
    plt.plot(CY_out, CY_out_pred, markerfacecolor='none', color = 'k'
             , marker='o', ls='', label='Data')
    plt.plot([min(CY_out),max(CY_out)],[min(CY_out),max(CY_out)]
             , color = 'r', ls = '-'
             , label = 'Parity Line')
    plt.xlabel("$C_{Y,out,meas}$ (M)")
    plt.xlim(left=0)
    plt.ylabel("$C_{Y,out,pred}$ (M)")
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.savefig('./png/example_11_5_2_parity.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_2_parity.png', dpi=300)
    plt.savefig('./pdf/example_11_5_2_parity.pdf')
    plt.show(block=False)

    # generate show and save residuals plots
    plt.figure(2) 
    plt.plot(Vdot*1.0E3, epsilon_expt, markerfacecolor='none', color = 'k'
             , marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\dot{V}$ (cm^3^ min^-1^)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_2_Vdot_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_2_Vdot_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_2_Vdot_residual.pdf')
    plt.show(block=False)

    plt.figure(3) 
    plt.plot(T_K, epsilon_expt, markerfacecolor='none', color = 'k'
             , marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("T (K)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_2_T_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_2_T_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_2_T_residual.pdf')
    plt.show(block=False)

    plt.figure(4) 
    plt.plot(CA_in, epsilon_expt, markerfacecolor='none', color = 'k'
             , marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{A,in}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_2_CA_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_2_CA_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_2_CA_residual.pdf')
    plt.show(block=False)

    plt.figure(5) 
    plt.plot(CB_in, epsilon_expt, markerfacecolor='none', color = 'k'
             , marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{B,in}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_2_CB_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_2_CB_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_2_CB_residual.pdf')
    plt.show(block=False)

    plt.figure(6) 
    plt.plot(CY_in, epsilon_expt, markerfacecolor='none', color = 'k'
             , marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{Y,in}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_2_CY_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_2_CY_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_2_CY_residual.pdf')
    plt.show(block=False)

    plt.figure(7) 
    plt.plot(CZ_in, epsilon_expt, markerfacecolor='none', color = 'k'
             , marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{Z,in}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_2_CZ_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_2_CZ_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_2_CZ_residual.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    