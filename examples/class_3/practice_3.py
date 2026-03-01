'''Calculations for the REB Class 3 Practice Assignment'''

import numpy as np
import pandas as pd
from reb_utils import Arrhenius_parameters
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi=300)

# given and known constants
R = 1.987 # BTU/lbmol/°R

# read experimental data and convert temperatures to °R
df = pd.read_csv('class_3_practice_data.csv')
df.columns = ['T','k']
T = df['T'].to_numpy() + 459.67 # °R
k_meas = df['k'].to_numpy() # /s

# calculate the Arrhenius parameters and statistics
k0, k0_ci, E, E_ci, r_squared = Arrhenius_parameters(k_meas,T,R)

# generate, show, and save an Arrhenius plot
k_AE = k0*np.exp(-E/R/T)
plt.figure()
plt.plot(1/T,np.log(k_meas),color='k',marker='o', ls='none'
         , label='Experimental Data')
plt.plot(1/T,np.log(k_AE),color='r', label='Arrhenius Expression')
plt.xlabel('T$^{-1}$')
plt.ylabel('ln k')
plt.legend()
plt.savefig('Python_Arrhenius_plot.png')
plt.savefig('Python_Arrhenius_plot.pdf')
plt.show()

# tabulate, display, and save the results
data = [['k0', f'{k0:.5g}', 's^-1^'],
    ['k0_lower_limit', f'{k0_ci[0]:.5g}', 's^-1^'],
    ['k0_upper_limit', f'{k0_ci[1]:.5g}', 's^-1^'],
    ['E', f'{E:.5g}', 'BTU lbmol^-1^'],
    ['E_lower_limit', f'{E_ci[0]:.5g}', 'BTU lbmol^-1^'],
    ['E_upper_limit', f'{E_ci[1]:.5g}', 'BTU lbmol^-1^^'],
    ['R_squared', f'{r_squared:.5g}', '']]
results = pd.DataFrame(data, columns=['item','value','units'])
print(' ')
print(results)
print(' ')
results.to_csv("class_3_practice_Python_results.csv", index=False)
