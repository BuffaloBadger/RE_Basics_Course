"""Calculations for the Class 13 Learning Activity in REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates

# global constants available to all function
# given
V = 4.430e3 # cm3
CS_in = 0.04 # g/cm3
CX_in = 0 # g/cm3
Vmax = 0.014 # /min
Km = 0.001 # g/cm3
Vdot_in = 60 # cm3/min

# cstr model function
def cstr_model_variables(init_guess):

    # solve the cstr design equations
    soln, success, message = solve_ates(cstr_residuals, init_guess)
    
    # check for solver issues
    if not success:
        print('')
        print(f"CSTR model function issue: {message}")
        print('')
        input('Press return to continue.')
    
    # return the cstr model variables
    return soln[0], soln[1] # mS, mX

# cstr residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    mS = guess[0]
    mX = guess[1]

    # calculate the additional unknowns
    Vdot = Vdot_in
    CS = mS/Vdot
    CX = mX/Vdot
    r1 = Vmax*CS*CX/(Km + CS)
    mS_in = CS_in*Vdot_in
    mX_in = CX_in*Vdot_in

    # evaluate the derivatives
    epsilon_1 = mS_in - mS - 2.2*V*r1
    epsilon_2 = mX_in - mX + V*r1

    # return the derivatives
    return np.array([epsilon_1, epsilon_2])

# deliverables function
def deliverables():
    # define an initial guess for the cstr model variables
    init_guess = np.array([0.5*CS_in, 0.1*CS_in]) *Vdot_in

    # solve the cstr design equations for the first inlet flow rate
    mS, mX = cstr_model_variables(init_guess)

    # tabulate, show and save the results
    data = [["CS", mS/Vdot_in, "g/cm3"]
            ,["CX", mX/Vdot_in, "g/cm3"]]
    results_df = pd.DataFrame(data, columns=["item","value","units"])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_13_ss_results.csv',index=False)

    # check the solution
    guess = np.array([mS, mX])
    residuals = cstr_residuals(guess)
    print('')
    print('Residuals Check:')
    print(residuals)
    print('')

# execution command
if __name__=="__main__":
    deliverables()
