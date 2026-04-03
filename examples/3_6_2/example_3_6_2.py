"""Calculations for REB Chapter 3, Example 2"""

# import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# set the dpi for figures
plt.rc("savefig", dpi=300)

# given and known constants
k1 = 3.37 # lbmol/h/ft^3/atm^-0.55
K1 = 12.01
f_CO_eq = 96.28
y_H2O_0 = 0.75
y_CO_0 = 0.24
y_CO2_0 = 0.01
P = 10 # atm
T = 675 # K
R = 0.08206 # L atm/mol/K

# basis
V = 1.0 # L

# initial moles, eqns 4 - 7
n_CO_0 = y_CO_0*P*V/R/T
n_H2O_0 = y_H2O_0*P*V/R/T
n_CO2_0 = y_CO2_0*P*V/R/T
n_H2_0 = 0

# create vectors for plot data
f_CO = np.linspace(0.0, 1.0, 100)
r_m0 = np.zeros(100, dtype=float, order='C')
r_m1 = np.zeros(100, dtype=float, order='C')

# calculate the predicted rates at each conversion
for i in range(0, np.size(f_CO)):
    # calculated the moles of CO, eqn 8
    n_CO = n_CO_0*(1-f_CO[i])
    # calculate the apparent extent, eqn 9
    extent = n_CO_0 - n_CO
    # calculate molar amounts, eqns 10 - 12
    n_H2O = n_H2O_0 - extent
    n_CO2 = n_CO2_0 + extent
    n_H2 = n_H2_0 + extent
    # calculate the partial pressures, eqns 13 - 16
    P_CO = n_CO*R*T/V
    P_H2O = n_H2O*R*T/V
    P_CO2 = n_CO2*R*T/V
    P_H2 = n_H2*R*T/V
    # calculate the rates
    r_m0[i] = k1*P_CO**0.9*P_H2O**0.25*P_CO2**-0.6
    factor = 1 - P_CO2*P_H2/K1/P_CO/P_H2O
    r_m1[i] = factor * r_m0[i]

# plot the results
plt.figure(1) # with conversions between 0 and 100%
plt.plot(100*f_CO,r_m0, color = 'r', label = 'without equilibrium factor')
plt.plot(100*f_CO,r_m1, color = 'b', label = 'with equilibrium factor')
plt.axhline(y=0, color = 'k')
plt.xlabel("CO Conversion (%)")
plt.ylabel("Predicted Rate (lbmol h$^{-1}$ ft$^{-3}$)")
plt.legend()
# save and show the figure
plt.savefig('example_3_6_2_fig_1.pdf')
plt.savefig('example_3_6_2_fig_1.png')
plt.show()

plt.figure(2) # with conversions between 85 and 100%
plt.plot(100*f_CO[84:],r_m0[84:], color = 'r', label = 'without equilibrium factor')
plt.plot(100*f_CO[84:],r_m1[84:], color = 'b', label = 'with equilibrium factor')
plt.axhline(y=0, color = 'k')
plt.axvline(f_CO_eq, color = 'k', ls = ':')
plt.xlabel("CO Conversion (%)")
plt.ylabel("Predicted Rate (lbmol h$^{-1}$ ft$^{-3}$)")
plt.legend()
# save and show the figure
plt.savefig('example_3_6_2_fig_2.pdf')
plt.savefig('example_3_6_2_fig_2.png')
plt.show()

# save the plot data to a .csv file
plot_data = pd.DataFrame({'Conversion':100*f_CO[0:98], 'r2':r_m0[0:98], 'r3':r_m1[0:98]})
plot_data.to_csv("example_3_6_2_results.csv", index=False)