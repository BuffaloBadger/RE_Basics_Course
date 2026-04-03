"""Calculations for REB Example 6.6.1"""

import example_6_6_1
import numpy as np
import pandas as pd

# read in the results from the assignment
results_df = pd.read_csv('example_6_6_1_results.csv')
# drop all but the first 3 rows
results_df = results_df.iloc[[0,1,2],:]

# set the initial guess
init_guess = np.array([example_6_6_1.nA_in, example_6_6_1.nB_in, 0.0, 0.0
            , example_6_6_1.T_in + 10.0])

# solve the design equations for a lower inlet temperature
example_6_6_1.T_in = 325 # K
solution = example_6_6_1.reactor_model_variables(init_guess)
# extract individual results
nA = solution[0]
nD = solution[2]
nU = solution[3]
T = solution[4]

# calculate the other quantities of interest
fA = 100*(example_6_6_1.nA_in - nA)/example_6_6_1.nA_in
S_D_U = nD/nU

# add the new results
results_df.loc[3] = ['Conversion(325 K)',f'{fA}','%']
results_df.loc[4] = ['Selectivity(325 K)',f'{S_D_U}','mol D per mol U']
results_df.loc[5] = ['T out (325)',f'{T}', 'K']

# solve the design equations for a higher inlet temperature
example_6_6_1.T_in = 375 # K
solution = example_6_6_1.reactor_model_variables(init_guess)
# extract individual unknowns
nA = solution[0]
nD = solution[2]
nU = solution[3]
T = solution[4]

# calculate the other quantities of interest
fA = 100*(example_6_6_1.nA_in - nA)/example_6_6_1.nA_in
S_D_U = nD/nU

# add the new results
results_df.loc[6] = ['Conversion(375 K)',f'{fA}','%']
results_df.loc[7] = ['Selectivity(375 K)',f'{S_D_U}','mol D per mol U']
results_df.loc[8] = ['T out (375)',f'{T}', 'K']


# display the results
print(results_df)

# save the results
results_df.to_csv('example_6_6_1_results.csv', index=False)
