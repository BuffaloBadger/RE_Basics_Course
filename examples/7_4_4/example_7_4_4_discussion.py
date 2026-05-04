"""Calculations for discussion of Example 7.4.4 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import example_7_4_4

# set graph resolution
plt.rc('savefig', dpi=300)

# read and set the optimum coolant flow rate
results_df = pd.read_csv('example_7_4_4_results.csv')
example_7_4_4.g_Vdot_ex=results_df.iat[0,1]

print("")
print(f"coolant flow rate: {results_df.iat[0,1]}")
print('')

# solve the reactor design equations using the optimum coolant flow
t_1, nA_1, nZ_1, T_1, Tex_1 = example_7_4_4.stage1_bstr_model_variables()
t_2, nA_2, _, T_2, _ = example_7_4_4.stage2_bstr_model_variables(t_1[-1], nA_1[-1], nZ_1[-1]
                                            , T_1[-1], Tex_1[-1])

# combine the profiles
t = np.concatenate((t_1, t_2))
nA = np.concatenate((nA_1, nA_2))
T = np.concatenate((T_1, T_2))

# calculate the conversion and instantaneous rate
r_inst = example_7_4_4.k0_1*np.exp(-example_7_4_4.E_1/example_7_4_4.R/T)\
    *nA/example_7_4_4.V*1000

# calculate and plot instantaneous rate profile for discussion
plt.figure(1) # instantaneous rate profile
plt.plot(t,r_inst)
plt.xlabel("Reaction Time (min)")
plt.ylabel("Instantaneous Rate (mmol/L/min)")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig('example_7_4_4_inst_rate_profile.png')
plt.savefig('example_7_4_4_inst_rate_profile.pdf')
plt.savefig('../../../RE_Basics/solutions/ch7_ex4/example_7_4_4_inst_rate_profile.png')
plt.show()