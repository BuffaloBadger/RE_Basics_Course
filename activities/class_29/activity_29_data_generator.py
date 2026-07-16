'''Data Generation for the Class 29 Learning Activity from REB, The Course'''

import numpy as np
import pandas as pd
import random
from reb_utils import solve_ivodes
from great_tables import GT, html
import os

# Constant inputs
Rpv = 82.06 # cc atm/mol/K
V = 100 # cc
P0 = 6.0 # atm

# kinetics parameters
k = 1.62E-6 # mol/cc/atm^1.5/min
T_ref = 275.0 + 273.15
E = 14000 # cal/mol
k0 = k*np.exp(E/1.987/T_ref)

# Define rate expression
def rate(P_A, P_B, temp):
    if P_A <= 0.0:
        reaction_rate = 0.0
    elif P_B <= 0.0:
        reaction_rate = 0.0
    else:
        reaction_rate = k0*np.exp(-E/1.987/temp)*P_A*np.sqrt(P_B)
    return reaction_rate

# Adjusted inputs
T_expt = np.array([225.0 + 273.15, 250 + 273.15, 275 + 273.15])
PAin = np.array([2.0, 3.0, 4.0])
t_reaction = np.arange(1.0, 25, 1.0)

# Create empty dataframe for the results
df = pd.DataFrame(columns=["T", "PA0", "t_meas", "P_meas"])

# Calculate the responses
for T in T_expt:
    for PA0 in PAin:
        PB0 = P0 - PA0
        for time in t_reaction:
            # Define the mole balances
            def mole_balances(t,n):
                PA = n[0]*Rpv*T/V
                PB = n[1]*Rpv*T/V
                r = rate(PA,PB,T)
                ddt = np.array([-r*V, -r*V, r*V])
                return ddt
            
            # Solve the mole balances
            t0 = 0
            n0 = [PA0*V/Rpv/T, PB0*V/Rpv/T, 0.0]
            f_var = 0
            f_val = time
            t, dep, success, message = solve_ivodes(t0, n0, f_var, f_val, mole_balances, False)

            # calculate the response
            nA = dep[0,-1]
            nB = dep[1,-1]
            nZ = dep[2,-1]
            P = (nA + nB + nZ)*Rpv*T/V

            # add +/- 0.05 atm random "error"
            random_error = (2*random.random() - 1.0)*0.05
            P = P + random_error

            # round to 2 decimal places
            P = round(P,2)

            # append the result to the dataframe
            df.loc[len(df.index)] = [T - 273.15, PA0, time, P]

# show and save the results
print('')
print(df)
df.to_csv('activity_29_data.csv',index=False)

# make sure a folder for saving png files exists
if not os.path.isdir('png'):
    os.makedirs('./png')

# generate a table showing the first 6 data points.
table_df = df.head(6)
data_tbl = (GT(table_df)
    .cols_label(T=html('T (°C)'),
                PA0=html('P<sub>A,0</sub> (atm)'),
                t_meas=html('t<sub>meas</sub> (min)'),
                P_meas=html('P<sub>meas</sub> (atm)')
    )
    .cols_align(
        align='center',
        columns=[0,1,2,3]
    )
)
data_tbl.show()
data_tbl.gtsave('png/data_table.png')

