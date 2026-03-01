'''Calculations for the REB Class 3 Learning Activity'''

# import libraries
from reb_utils import Arrhenius_parameters
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# set the dpi for figures
plt.rc("savefig", dpi=300)

# given and known constants
k01f = 6.97E6 # L/mol/min
E1f = 57.3 # kJ/mol
Ren = 8.314E-3 # kJ/mol/K
Rpv = 0.082057 # L atm/mol/K

# generate k_prime vs T data
T = np.linspace(273, 373, 100) # K
k_prime = k01f*np.exp(-E1f/(Ren*T))/(Rpv*T)**2

# calculate the Arrhenius parameters and statistics
k0, k0_ci, E, E_ci, r_squared = Arrhenius_parameters(k_prime,T,Ren)

# generate, show, and save an Arrhenius plot
k_AE = k0*np.exp(-E/Ren/T) # k prime calculated using the estimated Arrhenius parameters
plt.figure()
plt.plot(1/T,np.log(k_prime),color='k',marker='o', ls='none'
         , label='Apparent k')
plt.plot(1/T,np.log(k_AE),color='r', label='Arrhenius Expression')
plt.xlabel('T$^{-1}$')
plt.ylabel('ln k')
plt.legend()
plt.savefig('activity_3_fig_1.png')
plt.savefig('activity_3_fig_1.pdf')
plt.show()

# tabulate, display, and save the results
data = [['k0', f'{k0:.5g}', 'mol L^-1^ atm^-2^ min^-1^'],
    ['k0_lower_limit', f'{k0_ci[0]:.5g}', 'mol L^-1^ atm^-2^ min^-1^'],
    ['k0_upper_limit', f'{k0_ci[1]:.5g}', 'mol L^-1^ atm^-2^ min^-1^'],
    ['E', f'{E:.5g}', 'kJ mol^-1^'],
    ['E_lower_limit', f'{E_ci[0]:.5g}', 'kJ mol^-1^'],
    ['E_upper_limit', f'{E_ci[1]:.5g}', 'kJ mol^-1^'],
    ['R_squared', f'{r_squared:.5g}', '']]
results = pd.DataFrame(data, columns=['item','value','units'])
print(' ')
print(results)
print(' ')
results.to_csv("activity_3_results.csv", index=False)
