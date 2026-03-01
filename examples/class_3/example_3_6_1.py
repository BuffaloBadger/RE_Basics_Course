'''Calculations for REB Chapter 3, Example 1'''

# import libraries
import pandas as pd

# given and known constants
rho_bed = 155.0 * 0.0283 # lbm/m^3
k0_NH3 = 1.54e15 # kmol NH3 / m^3 / h

# Calculate k0 for generalized rate normalized per pound
k0 = k0_NH3/2/rho_bed

# tabulate, display, and save the results
data = [['k0', f'{k0:.3g}', 'kmol NH3 /lbm /h']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv("example_3_6_1_results.csv", index=False)
