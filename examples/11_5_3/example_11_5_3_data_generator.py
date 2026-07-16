import numpy as np
import pandas as pd
import scipy as sp
import random
import os

# Constant inputs
V = 3.0 # gal
R = 1.987 # BTU/lbmol/°R

# kinetics parameters
T_base = 110 + 459.7
kf = 15.98 # gal/lbmol/min
Ef = 12000 # BTU/lbmol
k0f = kf*np.exp(Ef/R/T_base)
print(f'k0f: {k0f:.3g}')
kr = 5.02 # gal/lbmol/min
Er = 19000 # BTU/lbmol
k0r = kr*np.exp(Er/R/T_base)
print(f'k0r: {k0r:.3g}')

# Define rate expression
def rate(kf, kr, concA, concY, concZ):
    return kf*concA**2 - kr*concY*concZ

# Adjusted inputs
T_levels = np.array([90, 110, 130]) + 459.7 # °R
VFR = np.array([0.5, 1.0, 1.5, 2.0, 2.5]) # gal/min
CA0 = np.array([0.01, 0.02, 0.03]) # lbmol/gal
CY0 = np.array([0.0, 0.01, 0.02, 0.03]) # lbmol/gal
CZ0 = np.array([0.0, 0.01, 0.02, 0.03]) # lbmol/gal

# Create empty dataframe for the results
df = pd.DataFrame(columns=['T (F)', 'Vdot (gal/min)', 'CA_in (lbmol/gal)'
                  , 'CY_in (lbmol/gal)', 'CZ_in (lbmol/gal)', 'CA_out (lbmol/gal)'])

# Calculate the responses
for T in T_levels:
    kf = k0f*np.exp(-Ef/R/T)
    kr = k0r*np.exp(-Er/R/T)
    for Vdot in VFR:
        for CAin in CA0:
            for CYin in CY0:
                for CZin in CZ0:
                    # define the mole balance
                    def mole_balance(fA):
                        nAin = CAin*Vdot
                        nYin = CYin*Vdot
                        nZin = CZin*Vdot
                        nA = nAin*(1-fA)
                        nY = nYin + nAin*fA
                        nZ = nZin + nAin*fA
                        CA = nA/Vdot
                        CY = nY/Vdot
                        CZ = nZ/Vdot
                        r = rate(kf, kr, CA, CY, CZ)
                        return nAin - nA - V*r
                    
                    # solve the mole balance
                    guess = [0.1]
                    soln = sp.optimize.root(mole_balance, guess)
                    if not(soln.success):
                        print("A solution was NOT obtained:")
                        print(soln.message)
                    fA = soln.x[0]

                    # add +/- 0.01 random "error"
                    random_error = (2*random.random() - 1.0)*0.05
                    fA += random_error

                    CAout = CAin*(1-fA)

                    # round to 2 decimal places
                    CAout = round(CAout,4)

                    # append the result to the dataframe
                    df.loc[len(df.index)] = [T - 459.7,Vdot,CAin,CYin
                                             ,CZin,CAout]

# display the results
print('\n')
print(df)

# make sure a REB folder exists for saving the results
reb_book_path = '../../../RE_Basics/solutions/ch11_ex3/'
if not os.path.isdir(reb_book_path):
    # create the folder
    os.makedirs(reb_book_path)

# save the results
df.to_csv('example_11_5_3_data.csv',index=False)
df.to_csv(reb_book_path + 'example_11_5_3_data.csv',index=False)