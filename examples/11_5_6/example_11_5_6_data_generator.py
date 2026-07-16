"""Data Generation for Example 11.5.6 from REB, The Book"""

import numpy as np
import pandas as pd
import random
from reb_utils import solve_ivodes
import os
from great_tables import GT, md, html

# Constant inputs
V = 1.0 # L
R = 8.314E-3 # kJ/mol/K

# Define rate expression
def rate(conc_A, temp):
    k0 = 4.23E8
    E = 68.0
    reaction_rate = k0*np.exp(-E/R/temp)*conc_A
    return reaction_rate

# Adjusted inputs
T_expt = np.array([65, 73, 82, 90]) + 273.15
CA0_expt = np.array([0.5, 1.0, 1.5])
t_reaction = np.array([5.0, 10.0, 15.0, 20.0, 25.0, 30.0])

# Create empty dataframe for the results
df = pd.DataFrame(columns=["Experiment","T", "CA0", "t_meas", "CA_meas"])

# Calculate the responses
expt_number = 0
for CA0 in CA0_expt:
    for T in T_expt:
        expt_number += 1
        for time in t_reaction:
            # Define the mole balances
            def mole_balances(t,C):
                ddt = [-rate(C,T)]
                return ddt
            
            # Solve the mole balances
            t0 = 0
            f_var = 0
            t, dep, success, message = solve_ivodes(t0, [CA0], f_var, time
                    , mole_balances, odes_are_stiff=False)

            # calculate the response
            CA = dep[0,-1]

            # add +/- 0.01 M random "error"
            random_error = (2*random.random() - 1.0)*0.08
            CA = CA + random_error

            # round to 2 decimal places
            CA = round(CA,2)

            # append the result to the dataframe
            df.loc[len(df.index)] = [expt_number, T - 273.15, CA0, time, CA]

# show and save the full data file
print('')
print(df)
df.to_csv('example_11_5_6_data.csv',index=False)

# make a table with the first six rows
table_df = df.head(6)
data_tbl = (GT(table_df)
    .tab_header(
        title='Results from Isothermal BSTR Experiments',
        subtitle='The first 6 of 72 data points are shown'
    )
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
                CA_meas=html('C<sub>A,f</sub><br>(M)')
    )
    .cols_align(
        align='center',
        columns=[0,1,2,3,4]
    )
)

# make sure a folder for saving png files exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# save the data table
data_tbl.gtsave('png/data_table.png', zoom=4.0)

