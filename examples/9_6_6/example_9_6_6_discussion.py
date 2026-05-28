"""Calculations for Discussion of Example 9.6.6 from REB, The Book"""

# import libraries
import pandas as pd
import example_9_6_6

example_9_6_6.eta = 1.0
example_9_6_6.Pconv = 0.0

# solve the reactor design equations
z, nA, nZ, P = example_9_6_6.pfr_model_variables()

# calculate the other quantities of interest
L = z[-1]
P_out = P[-1]

# read in the results from the assignment
results_df = pd.read_csv('example_9_6_6_results.csv')

# add the new results
n_rows = len(results_df.index)
results_df.loc[n_rows] = ['effectiveness factor',1.0,'']
results_df.loc[n_rows+1] = ['reactor length',L,'cm']

# display the results
print(' ')
print(results_df)
print('')

# save the results
results_df.to_csv('example_9_6_6_results.csv', index=False)
    