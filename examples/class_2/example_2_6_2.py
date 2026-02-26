"""Calculations for Example 2.6.2 in REB, The Book"""
# import libraries
import pandas as pd

# given and known constants
y_ethane = 0.9
y_air = 0.1
frac_O2_air = 0.21
frac_N2_air = 0.79
f_O2 = 0.5
s_ethene_to_CO2 = 3
n_total = 1 # mol basis

# calculate the initial moles of each component
n_ethane0 = y_ethane * n_total
n_air0 = y_air * n_total
n_O20 = frac_O2_air * n_air0
n_N20 = frac_N2_air * n_air0
n_total0 = n_ethane0 + n_air0

# calculate the extent of each reaction
extent2 = 2*n_O20*f_O2/(2*s_ethene_to_CO2 + 7)
extent1 = n_O20*f_O2 - 3*extent2

# calculate final mole fracton of CO2
y_CO2 = (2*extent2) / (n_total0 + extent1)

# tabulate the results
data = [["Initial n Ethane", n_ethane0, "mol"],
        ["Initial n Air", n_air0, "mol"],
        ["Initial n O2", n_O20, "mol"],
        ["Initial n N2", n_N20, "mol"],
        ["Extent of Reaction 1", extent1, "mol"],
        ["Extent of Reaction 2", extent2, "mol"],
        ["Final y CO2", y_CO2, "-"]]
resultsTable = pd.DataFrame(data, columns=["Item", "Value", "Units"])

# display and save the results
print(' ')
print(resultsTable)
print(' ')
resultsTable.to_csv("example_2_6_2_results.csv", index=False)