"""Calculations for Example 11.5.5 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes
from reb_utils import fit_to_SR_data
import os

# global constants available to all functions
# given
P = 1.0 # atm
m_cat = 3.0 # g
nTot_in = 2.5 / 60 # mol /min
K0 = 0.0132
dH = -9096 # cal /mol
# known
R = 1.987 # cal/mol/K

# experimental data as arrays
expt_df = pd.read_csv('example_11_5_5_data.csv')
yA_in = expt_df['yA_in'].to_numpy()
yB_in = expt_df['yB_in'].to_numpy()
yY_in = expt_df['yY_in'].to_numpy()
yZ_in = expt_df['yZ_in'].to_numpy()
T = expt_df['T'].to_numpy() + 273.15
PA_out = expt_df['PA_out'].to_numpy()

# global variables for the current values of the parameters
global g_k0, g_E, g_T, g_a, g_b, g_y
g_k0 = float('nan')
g_E = float('nan')
g_T = float('nan')
g_a = float('nan')
g_b = float('nan')
g_y = float('nan')

# PFR model function
def pfr_model_variables(yA_in, yB_in, yY_in, yZ_in, T, k0, E, alphaA, alphaB
        , alphaY):
    
    # make the rate expression parameters available to the derivatives function
    global g_k0, g_E, g_T, g_a, g_b, g_y
    g_k0 = k0
    g_E = E
    g_T = T
    g_a = alphaA
    g_b = alphaB
    g_y = alphaY

    #initial values and stopping criterion
    ind0 = 0
    nA0 = yA_in*nTot_in
    nB0 = yB_in*nTot_in
    nY0 = yY_in*nTot_in
    nZ0 = yZ_in*nTot_in
    dep0 =np.array([nA0, nB0, nY0, nZ0])
    f_var = 0
    f_val = m_cat

    # solve the PFR design equations
    m, dep, success, message = solve_ivodes(ind0, dep0, f_var, f_val
            , pfr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not(success):
        print('')
        print(f"PFR model solver issue: {message}")
        print('')
        input('Press return to continue or CTRL-C to exit.')

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]

    # return the PFR model variables
    return m, nA, nB, nY, nZ

# PFR derivatives function
def pfr_derivatives (m, dep):
    nA = dep[0]
    nB = dep[1]
    nY = dep[2]
    nZ = dep[3]
    ntot = nA + nB + nY + nZ
    PA = nA/ntot*P
    PB = nB/ntot*P
    PY = nY/ntot*P
    PZ = nZ/ntot*P
    k = g_k0*np.exp(-g_E/(R*g_T))
    K = K0*np.exp(-dH/(R*g_T))
    r = k * (PA**g_a) * (PB**g_b) * (PY**g_y) * (1 - PY*PZ/(K*PA*PB))
    return np.array([-r, -r, r, r])

# predicted responses function
def predicted_responses(adj_inputs, beta, E, alphaA, alphaB, alphaY):
    # allocate storage for the responses
    PA_out_pred = np.ones_like(yA_in)*float('nan')

    # extract k0
    k0 = 10**beta

    # loop through the data points
    for iExpt, input in enumerate(adj_inputs):
        # solve the PFR design equations
        m, nA, nB, nY, nZ = pfr_model_variables(yA_in[iExpt], yB_in[iExpt]
            , yY_in[iExpt], yZ_in[iExpt], T[iExpt], k0, E, alphaA, alphaB
            , alphaY)
        
        # calculate the response
        nA_out = nA[-1]
        nB_out = nB[-1]
        nY_out = nY[-1]
        nZ_out = nZ[-1]
        PA_out_pred[iExpt] = nA_out/(nA_out + nB_out + nY_out + nZ_out)*P
    
    # return the responses
    return PA_out_pred

# deliverables function
def deliverables():
    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([yA_in, yB_in, yY_in, yZ_in, T]))

    # define a guess for the parameters
    par_guess = [0.0, 10000.0, 1.0, 1.0, -1.0]

    # estimate the parameters
    param, param_ci, r_squared, PA_out_pred = fit_to_SR_data(par_guess
            , adj_inputs , PA_out, predicted_responses, use_rel_error=False)
    
    # extract the results
    k0 = 10**param[0]
    k0_CI = 10**param_ci[0,:]
    E = param[1]
    E_CI = param_ci[1,:]
    alphaA = param[2]
    alphaA_CI = param_ci[2,:]
    alphaB = param[3]
    alphaB_CI = param_ci[3,:]
    alphaY = param[4]
    alphaY_CI = param_ci[4,:]

    # calculate the experiment residuals
    epsilon_expt = PA_out - PA_out_pred

    # calculate the rate coefficient units
    order = format(alphaA + alphaB + alphaY,'.2f')
    kUnits = 'mol /g /min /atm^-' + order

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
    reb_book_path = '../../../RE_Basics/solutions/ch11_ex5/'
    if not os.path.isdir(reb_book_path):
        # create the folder
        os.makedirs(reb_book_path)

    # tabulate, show and save the fitting results
    results = [['k0', f'{k0:.3g}',f'[{k0_CI[0]:.3g}, {k0_CI[1]:.3g}]',f'{kUnits}'],
        ['E', f'{E:.0f}',f'[{E_CI[0]:.0f}, {E_CI[1]:.0f}]',' cal/mol'],
        ['alphaA', f'{alphaA:.2g}',f'[{alphaA_CI[0]:.2g}, {alphaA_CI[1]:.2g}]',''],
        ['alphaB', f'{alphaB:.2g}',f'[{alphaB_CI[0]:.2g}, {alphaB_CI[1]:.2g}]',''],
        ['alphaY', f'{alphaY:.2g}',f'[{alphaY_CI[0]:.2g}, {alphaY_CI[1]:.2g}]',''],
        ['R_squared', f'{r_squared:.3f}', '','']]
    results_df = pd.DataFrame(results, columns=['Parameter','Value','95% CI','Units'])
    print(" ")
    print(results_df)
    results_df.to_csv("./csv/example_11_5_5_results.csv", index=False)
    results_df.to_csv(reb_book_path + "example_11_5_5_results.csv", index=False)

    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(PA_out, PA_out_pred, color = 'k', marker='o', ls=''
             , label = 'Data')
    plt.plot([min(PA_out),max(PA_out)],[min(PA_out),max(PA_out)]
             , color = 'r', ls = '-', label = 'Parity Line')
    plt.xlabel("$P_{A,out}$ (atm)")
    plt.ylabel("$P_{A,out,pred}$ (atm)")
    plt.legend
    plt.tight_layout()
    plt.savefig('./png/example_11_5_5_parity.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_5_parity.png', dpi=300)
    plt.savefig('./pdf/example_11_5_5_parity.pdf')
    plt.show(block=False)

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(yA_in, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$y_{A,in}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_5_A_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_5_A_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_5_A_residual.pdf')
    plt.show(block=False)

    plt.figure() 
    plt.plot(yB_in, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$y_{B,in}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_5_B_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_5_B_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_5_B_residual.pdf')
    plt.show(block=False)

    plt.figure() 
    plt.plot(yY_in, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$y_{Y,in}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_5_Y_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_5_Y_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_5_Y_residual.pdf')
    plt.show(block=False)

    plt.figure() 
    plt.plot(yZ_in, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$y_{Z,in}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_5_Z_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_5_Z_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_5_Z_residual.pdf')
    plt.show(block=False)

    plt.figure() 
    plt.plot(T-273.15, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("T (°C)")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('./png/example_11_5_5_T_residual.png', dpi=300)
    plt.savefig(reb_book_path + 'example_11_5_5_T_residual.png', dpi=300)
    plt.savefig('./pdf/example_11_5_5_T_residual.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    