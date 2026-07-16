import numpy as np
import pandas as pd
import scipy as sp
import random
from reb_utils import solve_ivodes
import os

# Constant inputs
V = 0.5 # L
R = 1.987E-3 # kcal/mol/K
Rpv = 0.08206 # L atm/mol/K

# kinetics parameters
k = 0.00176 # mol/L/atm^2/min
E = 21.7 # kcal/mol
T_avg = 500 + 273.15
k0 = k*np.exp(E/R/T_avg)
print(k0)

# Define rate expression
def rate(P_A, P_B, temp):
    reaction_rate = k0*np.exp(-E/R/temp)*P_A*P_B
    return reaction_rate

# Adjusted inputs
T_expt = np.array([475, 500, 525]) + 273.15
PA0_expt = np.array([0.5, 1.0, 1.5])
PB0_expt = np.array([0.5, 1.0, 1.5])
t_meas = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 5.0, 10.0])

# Create empty dataframe for the results
df = pd.DataFrame(columns=["T", "PA0", "PB0", "t_meas", "fA"])

# Calculate the responses
for T in T_expt:
    for PA0 in PA0_expt:
        nA0 = PA0*V/Rpv/T
        for PB0 in PB0_expt:
            nB0 = PB0*V/Rpv/T
            for time in t_meas:
                # Define the mole balances
                def mole_balances(t,n):
                    PA = n[0]*Rpv*T/V
                    PB = n[1]*Rpv*T/V
                    r = rate(PA, PB, T)
                    ddt = np.array([-r*V, -r*V])
                    return ddt
                
                # Solve the mole balances
                t0 = 0
                n0 = np.array([nA0, nB0])
                f_var = 0
                f_val = time
                t, dep, success, message = solve_ivodes(t0, n0, f_var, f_val, mole_balances, False)

                # calculate the response
                nA = dep[0,-1]
                fA = (PA0*V/Rpv/T - nA)/(PA0*V/Rpv/T)

                # add +/- 0.01 random "error"
                random_error = (2*random.random() - 1.0)*0.01
                fA = fA + random_error

                # round to 3 decimal places
                fA = round(fA,3)

                # append the result to the dataframe
                df.loc[len(df.index)] = [T - 273.15, PA0, PB0, time, fA]

# display the results
print("\n")
print(df)

# make sure a REB folder exists for storing results
reb_book_path = '../../../RE_Basics/solutions/ch11_ex4/'
if not os.path.isdir(reb_book_path):
    # create the folder
        os.makedirs(reb_book_path)


# save the results
filename = "example_11_5_4_data.csv"
df.to_csv(filename,index=False)
df.to_csv(reb_book_path + filename,index=False)
