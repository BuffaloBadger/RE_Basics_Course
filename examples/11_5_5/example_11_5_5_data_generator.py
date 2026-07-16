import numpy as np
import pandas as pd
import scipy as sp
import random
from reb_utils import solve_ivodes
import os

# Constant inputs
R = 1.987 # cal /mol /K
Vhat = 22.4 # L/mol
mCat = 3.0 # g
nTot_in = 2.5 # mol /h
nTot_in = nTot_in / 60. # mol /min
VFR_in = nTot_in * Vhat * 1000 # sccm
P = 1.0 # atm
E = 27400 # cal /mol
k0 = 0.36 * np.exp(E/(R*(500 + 273.15))) # mol /g /h
k0 = k0 /60 # mol /g /min
K0 = 0.0132
dH = -9096 # cal /mol

# kinetics parameters
alpha_A = 0.9
alpha_B = 0.25
alpha_Y = -0.6
print('')
print(f'data generated with k0 = {k0:.4g} mol /g /min')
print(f'  E = {E:.4g} cal /mol')
print(f'  alpha_A = {alpha_A:.2f}')
print(f'  alpha_B = {alpha_B:.2f}')
print(f'  alpha_Y = {alpha_Y:.2f}')
print(f'VFR_in = {VFR_in:.3g} sccm')

# Define rate expression
def rate(P_A, P_B, P_Y, P_Z, T):
    K = K0*np.exp(-dH/R/T)
    k = k0*np.exp(-E/R/T)
    reaction_rate = k * (P_A**alpha_A) * (P_B**alpha_B) \
        * (1 - P_Y*P_Z/(K*P_A*P_B))
    if (P_Y > 0):
        reaction_rate = reaction_rate * (P_Y**alpha_Y)
    return reaction_rate

def derivatives(m,n):
    n_tot = np.sum(n)
    PA = n[0]/n_tot
    PB = n[1]/n_tot
    PY = n[2]/n_tot
    PZ = n[3]/n_tot
    r = rate(PA, PB, PY, PZ, T)
    return np.array([-r, -r , r, r])

# Adjusted inputs
y1 = np.array([0.3, 0.4, 0.5, 0.6, 0.7])
y2 = np.array([0.33, 0.5, 0.66])
T_range = np.array([350, 400, 450]) + 273.15

# integration constants
ind_0 = 0
f_var = 0
f_val = mCat

# Create empty dataframe for the results
df = pd.DataFrame(columns=['yA_in','yB_in','yY_in','yZ_in','T','PA_out'])

# initialize counter
count = 0
'''
# lowest values case
T = 350 + 273.15
yA_in = 0.3
yB_in = 0.33*(1 - yA_in)
yY_in = 0.33*(1 - yA_in - yB_in)
yZ_in = 1 - yA_in - yB_in - yY_in

# solve the design equations
ind_0 = 0
dep_0 = np.array([yA_in, yB_in, yY_in, yZ_in]) * nTot_in

f_var = 0
f_val = mCat

ind, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
        , derivatives, odes_are_stiff=False)

nA_out = dep[0,-1]
nB_out = dep[1,-1]
nY_out = dep[2,-1]
nZ_out = dep[3,-1]
nTot_out = nA_out + nB_out + nY_out + nZ_out
PA_out = nA_out/nTot_out*P
PB_out = nB_out/nTot_out*P
PY_out = nY_out/nTot_out*P
PZ_out = nZ_out/nTot_out*P

print('')
print(f'low value case: PA = {PA_out:.3g}, PB = {PB_out:.3g}, PY = {PY_out:.3g}, PZ = {PZ_out:.3g}')

# highest values case
yA_in = 0.7
yB_in = 0.66*(1 - yA_in)
yY_in = 0.66*(1 - yA_in - yB_in)
yZ_in = 1 - yA_in - yB_in - yY_in

dep_0 = np.array([yA_in, yB_in, yY_in, yZ_in]) * nTot_in
ind, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
        , derivatives, odes_are_stiff=False)

nA_out = dep[0,-1]
nB_out = dep[1,-1]
nY_out = dep[2,-1]
nZ_out = dep[3,-1]
nTot_out = nA_out + nB_out + nY_out + nZ_out
PA_out = nA_out/nTot_out*P
PB_out = nB_out/nTot_out*P
PY_out = nY_out/nTot_out*P
PZ_out = nZ_out/nTot_out*P

print('')
print(f'high value case: PA = {PA_out:.3g}, PB = {PB_out:.3g}, PY = {PY_out:.3g}, PZ = {PZ_out:.3g}')

'''

# Calculate the responses
for T in T_range:
    for yA_in in y1:
        for yB_in in y2:
            yB_in = yB_in*(1-yA_in)
            for yc in y2:
                yY_in = yc*(1-yA_in-yB_in)
                yZ_in = 1 - yA_in - yB_in - yY_in

                K = K0*np.exp(-dH/R/T)
                equil_factor = yY_in*yZ_in/K/yA_in/yB_in
                if equil_factor > 1.0:
                    count += 1
                
                # set the initial values
                dep_0 = np.array([yA_in, yB_in, yY_in, yZ_in]) * nTot_in

                # solve the design equations
                ind, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                        , derivatives, odes_are_stiff=True)

                # calculate the response
                nA_out = dep[0,-1]
                nB_out = dep[1,-1]
                nY_out = dep[2,-1]
                nZ_out = dep[3,-1]
                nTot_out = nA_out + nB_out + nY_out + nZ_out
                PA_out = nA_out/nTot_out

                # add +/- 0.01 random "error"
                random_error = (2*random.random() - 1.0)*0.01
                if PA_out + random_error > yA_in:
                    PA_out = yA_in
                else:
                    PA_out = PA_out + random_error

                # round to 3 decimal places
                yA_in = round(yA_in,3)
                yB_in = round(yB_in,3)
                yY_in = round(yY_in,3)
                yZ_in = round(yZ_in,3)
                PA_out = round(PA_out,3)

                # append the result to the dataframe
                df.loc[len(df.index)] = [yA_in, yB_in, yY_in, yZ_in, T - 273.15
                        , PA_out]

print('')
print('There were %i experiments that started beyond equilibrium' %(count))
print('')

# check that A decreased
PA_in = df['yA_in'].to_numpy()
PA_out = df['PA_out'].to_numpy()
for i, out in enumerate(PA_out):
    if out > PA_in[i]:
        print(f'PA increased in experiment {i}')

# make sure a REB folder exists for saving the results
reb_book_path = '../../../RE_Basics/solutions/ch11_ex5/'
if not os.path.isdir(reb_book_path):
    # create the folder
    os.makedirs(reb_book_path)


# show and save the results
print('')
print(df)
df.to_csv('example_11_5_5_data.csv',index=False)
df.to_csv(reb_book_path + 'example_11_5_5_data.csv',index=False)
