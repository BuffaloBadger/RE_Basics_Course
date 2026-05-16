"""Calculations for discussion of Example 7.4.3 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import example_7_4_3

# set graph resolution
plt.rc('savefig', dpi=300)

# read and set the optimum coolant flow rate
results_df = pd.read_csv('example_7_4_3_results.csv')

print("")
print(f"coolant flow rate: {results_df.iat[0,1]}")
print('')

# solve the reactor design equations using the optimum coolant flow
t, nA, nZ, T, Tex = example_7_4_3.bstr_model_variables(results_df.iat[0,1])

# calculate the conversion and instantaneous rate
r_inst = example_7_4_3.k0_1*np.exp(-example_7_4_3.E_1/example_7_4_3.R/T)\
    *nA/example_7_4_3.V*1000

# calculate and plot instantaneous rate profile for discussion
plt.figure(1) # instantaneous rate profile
plt.plot(t,r_inst)
plt.xlabel("Reaction Time (min)")
plt.ylabel("Instantaneous Rate (mmol/L/min)")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig('example_7_4_3_inst_rate_profile.png')
plt.savefig('example_7_4_3_inst_rate_profile.pdf')
plt.savefig('../../../RE_Basics/solutions/ch7_ex3/example_7_4_3_inst_rate_profile.png')
plt.show()