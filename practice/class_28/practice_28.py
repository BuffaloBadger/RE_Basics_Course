"""Calculations for the Class 28 Practice Assignment from REB, the Course"""

#import libraries
import pandas as pd
import numpy as np
import scipy as sp
from reb_utils import solve_ates
from reb_utils import fit_to_SR_data
import matplotlib.pyplot as plt
import os

# global constants available to all functions
# given
P = 3.0 # atm
# basis
V = 1.0 # L
# known
Rpv = 0.08206 # L*atm/mol/K
Ren = 8.314E-3 # kJ/mol/K

# global variables for quantities that cannot be passed to the 
# residuals function
gT = float('nan')
gtau = float('nan')
gyA = float('nan')
gk0 = float('nan')
gE = float('nan')
galphaA = float('nan')
galphaB = float('nan')

# experimental data as arrays
expt_df = pd.read_csv('practice_28_data.csv')
T = expt_df['T'].to_numpy() + 273.15 # K
tau = expt_df['tau'].to_numpy()
yAin = expt_df['yA_in'].to_numpy()
CZ_meas = (expt_df['CZ_out'].to_numpy())/1000. # mol/L

# CSTR model function
def cstr_model_variables(T, tau, yA, k0, E, alphaA, alphaB):
    # make adjusted inputs and rate expression parameters globally available
    global gT, gtau, gyA, gk0, gE, galphaA, galphaB
    gT = T
    gtau = tau
    gyA = yA
    gk0 = k0
    gE = E
    galphaA = alphaA
    galphaB = alphaB

    # guess the solution
    dotVin = V/tau
    nA_guess = yA*P*dotVin/(Rpv*T)
    nB_guess = (1 - yA)*P*dotVin/(Rpv*T)
    nZ_guess = 0.0
    initial_guess = [nA_guess, nB_guess, nZ_guess]
     
	# solve the ATEs
    soln, success, message = solve_ates(cstr_residuals, initial_guess)

    # check that the solution is converged
    if not(success):
        print('')
        print(f'CSTR model function issue: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit')

    # extract the unknowns
    nA_out = soln[0]
    nB_out = soln[1]
    nZ_out = soln[2]

    return nA_out, nB_out, nZ_out

# CSTR residuals function
def cstr_residuals(guess):
    # extract the indiviaual guesses
    nA_out = guess[0]
    nB_out = guess[1]
    nZ_out = guess[2]

    # calculate the additional unknowns
    Vdot_in = V/gtau
    nA_in = gyA*P*Vdot_in/(Rpv*gT)
    nB_in = (1-gyA)*P*Vdot_in/(Rpv*gT)
    PA = nA_out*P/(nA_out + nB_out + nZ_out)
    PB = nB_out*P/(nA_out + nB_out + nZ_out)
    k = gk0*np.exp(-gE/(Ren*gT))
    r = k*(PA**galphaA)*(PB**galphaB)

    # evaluate the residuals
    epsilon_1 = nA_in - nA_out - V*r
    epsilon_2 = nB_in - nB_out - V*r
    epsilon_3 = -nZ_out + V*r

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3])

# predicted responses function
def predicted_responses(adj_inputs, beta, E, alphaA, alphaB):
    # allocate storage for the responses
    CZ_pred = np.ones_like(CZ_meas)*float('nan')

    # extract the pre-exponential factor
    k0 = 10**beta

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # solve the reactor design equations
        T_K = input[0]
        tau = input[1]
        yA = input[2]
        nA_out, nB_out, nZ_out = cstr_model_variables(T_K, tau, yA, k0, E
            , alphaA, alphaB)
        
        # calculate the model-predicted response
        CZ_pred[i] = nZ_out*P/((nA_out + nB_out + nZ_out)*Rpv*T_K)

    # return the responses
    return CZ_pred

# deliverables function
def deliverables():
    # guess the parameters
    #par_guess = [0.0, 40.0, 1.0, 1.0] # original guess
    par_guess = [3.0, 100.0, 1.0, 1.0]

    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([T, tau, yAin]))

    # estimate the parameters
    param, param_ci, r_squared, CZ_pred = fit_to_SR_data(par_guess
        , adj_inputs, CZ_meas, predicted_responses, use_rel_error=False)

    # extract the parameter estimates and their confidence intervals
    beta = param[0]
    beta_CI = param_ci[0,:]
    k0 = 10**beta
    k0_CI = 10**beta_CI 
    E = param[1]
    E_CI = param_ci[1,:]
    alphaA = param[2]
    alphaA_CI = param_ci[2,:]
    alphaB = param[3]
    alphaB_CI = param_ci[3,:]

    # generate show and save the results
    results = [
        ['k0', f'{k0:.3g}', f'{k0_CI[0]:.3g}', f'{k0_CI[1]:.3g}'
            , 'mol L^-1^ atm^-2^ s^-1^'],
        ['E', f'{E:.3g}', f'{E_CI[0]:.3g}', f'{E_CI[1]:.3g}', 'kJ mol^-1^'],
        ['alphaA', f'{alphaA:.3g}', f'{alphaA_CI[0]:.3g}', f'{alphaA_CI[1]:.3g}'
            , ''],
        ['alphaB', f'{alphaB:.3g}', f'{alphaB_CI[0]:.3g}', f'{alphaB_CI[1]:.3g}'
            , ''],
        ['R_squared', f'{r_squared:.3g}', '','','']]
    results_df = pd.DataFrame(results, columns=['Parameter','Value'
        ,'lower_limit', 'upper_limit', 'units'])
    results_df.to_csv("results.csv", index=False)
    print(' ')
    print(results_df)

    # calculate experiment residuals for residuals plots
    epsilon_expt = CZ_meas - CZ_pred

    # make sure a png folder exists
    if not os.path.isdir('png'):
        os.makedirs('./png')
        
    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(CZ_meas, CZ_pred, color = 'b', marker='x', ls=''
             , label = 'Data')
    plt.plot([min(CZ_meas),max(CZ_meas)],[min(CZ_meas),max(CZ_meas)]
             , color = 'k', ls = '-', label = 'Parity Line')
    plt.xlabel("Measured $C_{Z,out}$ (M)")
    plt.ylabel("Predicted $C_{Z,out}$ (M)")
    plt.legend
    plt.tight_layout()
    plt.savefig('png/parity.png', dpi=300)
    plt.show(block=False)

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(T-273.15, epsilon_expt, color = 'b', marker='x', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("T (°C)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('png/T_residual.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(tau, epsilon_expt, color = 'b', marker='x', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$\\tau$ (s)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('png/tau_residual.png')
    plt.show(block=False)

    plt.figure() 
    plt.plot(yAin, epsilon_expt, color = 'b', marker='x', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$y_{A,0}$")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('png/yA_residual.png')
    plt.show()

if __name__=="__main__":
    deliverables()
