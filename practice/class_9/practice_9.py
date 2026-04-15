"""Calculations for Practice Assignment 9 from REB, The Course."""

# import libraries
import numpy as np
import pandas as pd
from reb_utils import solve_ates

# global constants
# given
T_in = 150 + 273.15; # K
CA_in = 2.0; # mol/L
V = 500.0; # L
Vdot_in = 250.0; # L/h
T_ex = 180 + 273.15; # K
A = 2.0; # m2
U = 500E3; # cal/m2/h/K
k01 = 1.14E9; # L/mol/h
E1 = 16200.; # cal/mol
Cp = 1.17E3; # cal/L/K
dH1 = 18.2E3; # cal/mol
# known
R = 1.987; # cal/mol/K
# calculated
nA_in = CA_in*Vdot_in;

# cstr reactor model function
def cstr_model_variables(init_guess):
    # solve the cstr design equations
    soln, success, message = solve_ates(cstr_residuals, init_guess)

    # check for solver issues
    if not success:
        print(f".  CSTR model function error: {message}")

    # extract and return the results
    return soln[0], soln[1], soln[2]

# cstr residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nZ = guess[1]
    T = guess[2]

    # calculate the additional unknowns
    Vdot = Vdot_in
    k1 = k01*np.exp(-E1/(R*T))
    CA = nA/Vdot
    r1 = k1*CA**2
    Qdot = U*A*(T_ex - T)

    # evaluate the residuals
    epsilon_1 = nA_in - nA - r1*V
    epsilon_2 = -nZ + r1*V
    epsilon_3 = Qdot - Vdot_in*Cp*(T - T_in) - r1*V*dH1
    return np.array([epsilon_1, epsilon_2, epsilon_3])

# deliverables function
def deliverables():
    # define the initial guess for the cstr reactor variables
    init_guess = np.array([nA_in, 0.0, T_in - 5])

    # solve the cstr model equations
    nA, nZ, T = cstr_model_variables(init_guess)

    # calculate the conversion
    fA = 100*(nA_in - nA)/nA_in

    # tabulate, display and save the results
    data = [["fA", f"{fA:.2f}","%"],
            ["T", f"{T - 273.15:.2f}","°C"]]
    results_df = pd.DataFrame(data, columns=["Item", "Value", "Units"])
    print(" ")
    print(results_df)
    results_df.to_csv("practice_9_results.csv", index=False)
    return

# execution command
if __name__ == "__main__":
    deliverables()
