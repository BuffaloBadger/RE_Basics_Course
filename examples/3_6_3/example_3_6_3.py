'''Calculations for REEB Chapter 3 Example 3'''

# import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# import the REB Arrhenius parameter utility function
from reb_utils import Arrhenius_parameters

# set the dpi for figures
plt.rc("savefig", dpi=300)

# given and known constants
R = 8.314E-3 # kJ/mol/K

# read the data file and convert the temperature to K
df = pd.read_csv('example_3_6_3_data.csv')
df.columns = ['i','T','k']
T = df['T'].to_numpy() + 273.15 # K
k_meas = df['k'].to_numpy() # L/mol/min

# calculate the Arrhenius parameters and statistics
k0, k0_ci, E, E_ci, r_squared = Arrhenius_parameters(k_meas,T,R)

# generate, show, and save an Arrhenius plot
k_pred = k0*np.exp(-E/R/T)
plt.figure()
plt.plot(1/T,np.log(k_meas),color='k',marker='o', ls='none'
        , label='Experimental Data')
plt.plot(1/T,np.log(k_pred),color='r', label='Arrhenius Expression')
plt.xlabel('T$^{-1}$ (K$^{-1}$)')
plt.ylabel('ln(k)')
plt.legend()
plt.savefig('example_3_6_3_fig_1.png')
plt.savefig('example_3_6_3_fig_1.pdf')
plt.show()

# tabulate, display, and save the results
data = [['k0', f'{k0:.3g}', 'L mol^-1^ min^-1^'],
    ['k0_lower_limit', f'{k0_ci[0]:.3g}', 'L mol^-1^ min^-1^'],
    ['k0_upper_limit', f'{k0_ci[1]:.3g}', 'L mol^-1^ min^-1^'],
    ['E', f'{E:.3g}', 'kJ mol^-1^'],
    ['E_lower_limit', f'{E_ci[0]:.3g}', 'kJ mol^-1^'],
    ['E_upper_limit', f'{E_ci[1]:.3g}', 'kJ mol^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
results = pd.DataFrame(data, columns=['item','value','units'])
print(" ")
print(results)
print(" ")
results.to_csv("example_3_6_3_results.csv", index=False)
