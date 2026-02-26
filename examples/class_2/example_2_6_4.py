"""Calculations for Example 2.6.3 in REB, The Book"""
# import libraries
import pandas as pd
import numpy as np

# given and known constants
n_CH4_0 = 1 # mol
n_H2O_0 = 2.5 # mol
T = 920 + 273.15 # K
f_CH4 = 0.88
y_CO = 0.65
dH_298 = np.array([-17.89, -57.8, -26.41, 0, -94.05]) # kcal /mol
a = np.array([4.598, 7.701, 7.373, 6.483, 4.728]) # cal /mol /K
b = np.array([1.245E-2, 4.595E-4, -3.07E-3, 2.215E-3, 1.754E-2]) # cal /mol /K^2
c = np.array([2.86E-6, 5.521E-6, 6.662E-6, -3.298E-6, -1.388E-5]) # cal /mol /K^3
d = np.array([-2.703E-9, -0.859E-9, -3.037E-9, 1.826E-9, 4.097E-9]) # cal /mol /K^4

# calculate the heat of reaction 1 at 298 K
dH_1_298 = - dH_298[0] - dH_298[1] + dH_298[2] + 3*dH_298[3] # kcal / mol

# calculate the coefficients in the expression for delta H1
a1 = - a[0] - a[1] + a[2] + 3*a[3]
b1 = (- b[0] - b[1] + b[2] + 3*b[3])/2
c1 = (- c[0] - c[1] + c[2] + 3*c[3])/3
d1 = (- d[0] - d[1] + d[2] + 3*d[3])/4

# calculate delta H1 at T
dH_1_T = dH_1_298*1000 + a1*(T - 298) + b1*(T**2 - 298**2) + c1*(T**3 - 298**3) + d1*(T**4 - 298**4) # cal / mol
dH_1_T = dH_1_T / 1000 # kcal / mol

# calculate the heat of reaction 2 at 298 K
dH_2_298 = - dH_298[0] - dH_298[4] + 2* dH_298[2] + 2*dH_298[3] # kcal / mol

# calculate the coefficients in the expression for delta H2
a2 = - a[0] - a[4] + 2*a[2] + 2*a[3]
b2 = (- b[0] - b[4] + 2*b[2] + 2*b[3])/2
c2 = (- c[0] - c[4] + 2*c[2] + 2*c[3])/3
d2 = (- d[0] - d[4] + 2*d[2] + 2*d[3])/4

# calculate delta H2 at T
dH_2_T = dH_2_298*1000 + a2*(T - 298) + b2*(T**2 - 298**2) + c2*(T**3 - 298**3) + d2*(T**4 - 298**4) # cal / mol
dH_2_T = dH_2_T / 1000 # kcal / mol

# tabulate the results
results = [["ΔH₁ (298 K)", dH_1_298, "kcal/mol"],
           ["a1", a1, "cal/mol/K"],
           ["b1", b1, "cal/mol/K²"],
           ["c1", c1, "cal/mol/K³"],
           ["d1", d1, "cal/mol/K⁴"],
           ["ΔH₁ (T)", dH_1_T, "kcal/mol"],
           ["ΔH₂ (298 K)", dH_2_298, "kcal/mol"],
           ["a2", a2, "cal/mol/K"],
           ["b2", b2, "cal/mol/K²"],
           ["c2", c2, "cal/mol/K³"],
           ["d2", d2, "cal/mol/K⁴"],
           ["ΔH₂ (T)", dH_2_T, "kcal/mol"]]

# calculate the final moles of CH4 and CO
n_CH4 = n_CH4_0 * (1-f_CH4) # mol
n_CO = y_CO * n_CH4_0 # mol

# calculate the extents of reactions
extent2 = (y_CO - f_CH4) * n_CH4_0 # mol
extent1 = f_CH4*n_CH4_0 - extent2 # mol

# calculate the net heat absorbed
netQ = -extent1*dH_1_T - extent2*dH_2_T # kcal

# add to the results table
added_results = [["Extent of Reaction 1", extent1, "mol"],
                 ["Extent of Reaction 2", extent2, "mol"],
                 ["Net Heat Released", netQ, "kcal"]]
results.extend(added_results)

# calculate the final moles of all species
n_CH4_final = n_CH4_0 - extent1 - extent2
n_CO_final = extent1 + 2*extent2
n_H2_final = 3*extent1 + 2*extent2
n_H2O_final = n_H2O_0 - extent1
n_CO2_final = - extent2

# calculate the final percent compositions
n_total_final = n_CH4_final + n_CO_final + n_H2_final + n_H2O_final + n_CO2_final
pct_CH4 = (n_CH4_final / n_total_final) * 100
pct_CO = (n_CO_final / n_total_final) * 100
pct_H2 = (n_H2_final / n_total_final) * 100
pct_H2O = (n_H2O_final / n_total_final) * 100
pct_CO2 = (n_CO2_final / n_total_final) * 100

# add to the results table
added_results_2 = [["Final # CH4", n_CH4_final, "mol"],
                   ["Final # CO", n_CO_final, "mol"],
                   ["Final # H2", n_H2_final, "mol"],
                   ["Final # H2O", n_H2O_final, "mol"],
                   ["Final # CO2", n_CO2_final, "mol"],
                   ["Final % CH4", pct_CH4, "%"],
                   ["Final % CO", pct_CO, "%"],
                   ["Final % H2", pct_H2, "%"],
                   ["Final % H2O", pct_H2O, "%"],
                   ["Final % CO2", pct_CO2, "%"]]
results.extend(added_results_2)
results_table = pd.DataFrame(results, columns=["Item", "Value", "Units"])

# display and save the results
print(' ')
print(results_table)
print(' ')
results_table.to_csv('example_2_6_4_results.csv', index=False)
