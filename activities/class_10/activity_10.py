"""Calculations for the Class 10 Learning Activity in REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates

# global constants
# given
V = 4 # L
k01 = 2.59e9 # min-1
E1 = 16500 # cal/mol
dH1 = -22200 # cal/mol
Cp = 440 # cal/L/K
Vex = 0.5 #L
A = 0.6 # ft2
U = 5.65e3/60 # cal/ft2/min/K
rho_ex = 1 # kg/L
Cp_ex = 1000 # cal/kg/K
CA_in = 2 # mol/L
T_in = 60 + 273.15 # K
Tex_in = 20 + 273.15 # K
m_ex = 0.2 # kg/min
tau_values = np.array([50, 60]) # min
# known
R = 1.987 # cal/mol/K

# allocate global storage for the current value of tau
global g_tau
g_tau = float('NaN')

# CSTR reactor model function
def cstr_model_variables(tau):
    # make tau available to the CSTR residuals function
    global g_tau
    g_tau = tau

    # define an initial guess for the cstr reactor model variables
    Vdot_in = V/tau
    nA_in = CA_in*Vdot_in
    guess = np.array([nA_in, 0, T_in + 5, Tex_in + 5])

    # solve the cstr design equations
    soln, success, message = solve_ates(cstr_residuals, guess)

    # check for solver issues
    if not success:
        print('')
        print("CSTR model function error: " + message)
        print('')
        input("Press enter to continue.")
    
    # return the solution
    return soln[0], soln[1], soln[2], soln[3] # nA, nZ, T, and Tex

# CSTR residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nZ = guess[1]
    T = guess[2]
    Tex = guess[3]

    # calculate the additional unknowns
    Vdot_in = V/g_tau
    Vdot = Vdot_in
    nA_in = CA_in*Vdot_in
    nZ_in = 0
    k1 = k01 * np.exp(-E1 / (R * T))
    CA = nA / Vdot
    r1 = k1 * CA
    Qdot = U * A * (Tex - T)

    # evaluate the residuals
    epsilon_1 = nA_in - nA - r1 * V
    epsilon_2 = nZ_in - nZ + r1 * V
    epsilon_3 = Vdot_in*Cp*(T - T_in) + V*r1*dH1 - Qdot
    epsilon_4 = -m_ex*Cp_ex*(Tex - Tex_in) - Qdot

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4])

# deliverables function
def deliverables():
    # allocate storage for the conversions and temperatures corresponding to each space time
    fA = np.ones_like(tau_values) * float('NaN')
    T = np.ones_like(tau_values) * float('NaN')
    Tex = np.ones_like(tau_values) * float('NaN')

    # loop through the space times
    for i, tau in enumerate(tau_values):
        # solve the cstr design equations
        nA, nZ, T[i], Tex[i] = cstr_model_variables(tau)

        # calculate the conversion
        Vdot_in = V/tau
        nA_in = CA_in*Vdot_in
        fA[i] = 100*(nA_in - nA)/nA_in

    # tabulate, display, and save the results
    results_df = pd.DataFrame({'Space Time (min)': tau_values, 'Conversion (%)':fA, 'T (°C)': T - 273.15
                               , 'Coolant T (°C)': Tex - 273.15})
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_10_results.csv', index=False)

# execution command
if __name__ == "__main__":
    deliverables()
