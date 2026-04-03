"""Calculations for Example 2.6.1 in REB, The Book"""
# import libraries
import pandas as pd

# given and known constants
y_CO = 0.5
y_air = 0.5
nu_CO = -1
nu_O2 = -0.5
frac_O2_air = 0.21
nDotTotalIn = 1.0

# calculate inlet molar flow rates
nDot_CO_in = y_CO * nDotTotalIn
nDot_air_in = y_air * nDotTotalIn
nDot_O2_in = frac_O2_air * nDot_air_in

# calculate stoichiometric equivalences
stoichEqCO = nDot_CO_in / (-nu_CO)
stoichEqO2 = nDot_O2_in / (-nu_O2)

# determine limiting reactant
if stoichEqCO < stoichEqO2:
    limitingReactant = 'CO'
elif stoichEqO2 < stoichEqCO:
    limitingReactant = 'O~2~'
else:
    limitingReactant = 'Neither, both are in exact stoichiometric amounts'

# tabulate the results
data = [["nDot CO in", nDot_CO_in, "mol s^-1^"],
        ["nDot Air in", nDot_air_in, "mol s^-1^"],
        ["nDot O2 in", nDot_O2_in, "mol s^-1^"],
        ["CO equivalence", stoichEqCO, "mol s^-1^"],
        ["O2 equivalence", stoichEqO2, "mol s^-1^"],
        ["Limiting Reactant", limitingReactant, ""]]
resultsTable = pd.DataFrame(data, columns=["Item", "Value", "Units"])

# display and save the results
print(' ')
print(resultsTable)
print(' ')
resultsTable.to_csv("example_2_6_1_results.csv", index=False)
