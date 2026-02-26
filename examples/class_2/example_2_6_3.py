"""Calculations for Example 2.6.3 in REB, The Book"""
# import libraries
import pandas as pd

# given and known constants
n_N2O5_in = 0.5 # mol/min
n_N2_in = 0.5 # mol/min
n_total_in = n_N2O5_in + n_N2_in
y_N2_out = 0.3

# calculate the extent of reaction
extent = (n_N2_in - y_N2_out*n_total_in) / (3*y_N2_out)

# calculate the outlet N2O5 flow rate
n_N2O5 = n_N2O5_in - 2*extent

# calculate the conversion
f_N2O5 = 100*(n_N2O5_in - n_N2O5) / n_N2O5_in

# tabulate the results
data = [["Extent", extent, "mol/min"],
        ["Outlet N2O5", n_N2O5, "mol/min"],
        ["N2O5 Conversion", f_N2O5, "%"]]
results_table = pd.DataFrame(data, columns=["Item", "Value", "Units"])

# display and save the results
print(' ')
print(results_table)
print(' ')
results_table.to_csv('example_2_6_3_results.csv', index=False)