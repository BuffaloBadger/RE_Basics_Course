'''Calculations for REB Chapter 3 Example 4'''

# import libraries
from reb_utils import Arrhenius_parameters
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# set the dpi for figures
plt.rc("savefig", dpi=300)

# given and known constants
sigma = 38.5 * 1E-20 # m^2
mu = 1.063E-25 # kg
E = 184 # kJ/mol
k_B = 1.381E-23 # J/K
N_Av = 6.022E23 # 1/mol
R = 8.314E-3 # kJ/mol/K

# generate k_CT vs T data
T = np.linspace(300.0, 400.0, 100) # K
k_CT = N_Av*sigma*np.sqrt(2*k_B*T/np.pi/mu)*np.exp(-E/(R*T)) # m^3/mol/s

# calculate the Arrhenius parameters and statistics
k0, k0_ci, E, E_ci, r_squared = Arrhenius_parameters(k_CT,T,R)

# generate, show, and save a model plot
k_AE = k0*np.exp(-E/R/T)
plt.figure()
plt.plot(1/T,np.log(k_CT),color='k',marker='o', ls='none'
         , label='Collision Theory k')
plt.plot(1/T,np.log(k_AE),color='r', label='Arrhenius Expression')
plt.xlabel('T$^{-1}$ (K$^{-1}$)')
plt.ylabel('k (m$^3$ mol$^{-1}$ s$^{-1}$)')
plt.legend()
plt.savefig('example_3_6_4_fig_1.png')
plt.savefig('example_3_6_4_fig_1.pdf')
plt.show()

# tabulate, display, and save the results
data = [['k0', f'{k0:.3g}', 'm^3^ mol^-1^ s^-1^'],
    ['k0_lower_limit', f'{k0_ci[0]:.3g}', 'm^3^ mol^-1^ s^-1^'],
    ['k0_upper_limit', f'{k0_ci[1]:.3g}', 'm^3^ mol^-1^ s^-1^'],
    ['E', f'{E:.3g}', 'kJ mol^-1^'],
    ['E_lower_limit', f'{E_ci[0]:.3g}', 'kJ mol^-1^'],
    ['E_upper_limit', f'{E_ci[1]:.3g}', 'kJ mol^-1^'],
    ['R_squared', f'{r_squared:.3g}', '']]
results = pd.DataFrame(data, columns=['item','value','units'])
print(' ')
print(results)
print(' ')
results.to_csv("example_3_6_4_results.csv", index=False)
