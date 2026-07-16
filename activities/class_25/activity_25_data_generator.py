import numpy as np
import pandas as pd
import scipy as sp
from reb_utils import solve_ivodes
import random
from great_tables import GT, html
import os

# Constant inputs
V = 50.0 # ml

# kinetics parameters
Vm = 1.15E-4 # mmol/ml/min
Km = 0.0021 # mmol/ml

# Define rate expression
def rate(CS):
    reaction_rate = Vm*CS/(Km + CS)
    return reaction_rate

# Adjusted inputs
CS0 = np.array([15.0, 10.0, 5.0]) * 1E-3 # mmol/mL
t_reaction = np.linspace(5.0, 120.0, 24) # min

# Create empty dataframe for the results
df = pd.DataFrame(columns=["CS0", "t_meas", "CP_meas"])

# Calculate the responses
for CS_init in CS0:
    for time in t_reaction:
        # Define the mole balances
        def mole_balances(t,n): #n[0] = nS, n[1] = nP, n[2] = nH2O
            C_S = n[0]/V
            r = rate(C_S)
            ddt = np.array([-r*V, r*V, r*V])
            return ddt
                
        # Solve the mole balances
        t0 = 0
        n0 = np.array([CS_init*V, 0.0, 0.0])
        f_var = 0
        f_val = time
        t, dep, success, message = solve_ivodes(t0, n0, f_var, f_val
                , mole_balances, odes_are_stiff=False)

        # calculate the response
        CP = dep[1,-1]/V

        # add +/- 0.01 random "error"
        random_error = (2*random.random() - 1.0)*0.0004
        CP = CP + random_error

        # convert to mmol/L
        CS = 1000.0*CS_init
        CP = 1000.0*CP

        # round CP to 2 decimal places
        CP = round(CP,2)

        # append the result to the dataframe
        df.loc[len(df.index)] = [CS, time, CP]

# display the results
print('')
print(df)

# save the results
df.to_csv('activity_25_data.csv',index=False)

# make sure a folder for saving png files exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# generate a table showing the first 6 data points.
table_df = df.head(6)
data_tbl = (GT(table_df)
    .cols_label(CS0=html('C<sub>S,0</sub><br>(mmol L<sup>-1</sup>)'),
                t_meas=html('t<sub>meas</sub><br>(min)'),
                CP_meas=html('C<sub>P,meas</sub><br>(mmol L<sup>-1</sup>)')
    )
    .cols_align(
        align='center',
        columns=[0,1,2]
    )
)
data_tbl.show()
data_tbl.gtsave('png/data_table.png')
