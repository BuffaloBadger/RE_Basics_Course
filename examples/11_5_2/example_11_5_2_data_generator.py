import numpy as np
import pandas as pd
import scipy as sp
import random
import os

# Constant inputs
V = 0.1 # L

# kinetics parameters
k0 = 5.29E9 # mol/L/min
E = 12100 # cal/mol

# Define rate expression
def rate(concA, concB, concZ,T):
    return k0*np.exp(-E/1.987/T)*concA*concB/concZ**2

# Adjusted inputs
Texpt = np.array([300, 320, 340, 360])
VFR = np.array([50, 100]) # cc/min
VFR = VFR/1000.0 # L/min
CA0 = np.array([0.5, 1.0, 2.5, 5]) # M
CB0 = np.array([0.5, 1.0, 2.5, 5]) # M
CY0 = np.array([0.5, 1.0, 2.5, 5]) # M
CZ0 = np.array([0.5, 1.0, 2.5, 5]) # M

# Create empty dataframe for the results
df = pd.DataFrame(columns=['T', 'Vdot', 'CA_in', 'CB_in', 'CY_in', 'CZ_in'
            , 'CY_out'])

# Calculate the responses
for T in Texpt:
    for Vdot in VFR:
        for CAin in CA0:
            for CBin in CB0:
                for CYin in CY0:
                    for CZin in CZ0:
                        # define the mole balance
                        def mole_balance(CYout):
                            extent = CYout - CYin
                            CA = CAin - extent
                            CB = CBin - extent
                            CZ = CZin + extent
                            r = rate(CA, CB, CZ, T)
                            return (CAin - CA)*Vdot - V*r
                        
                        # solve the mole balance
                        guess = [CYin + 0.5]
                        soln = sp.optimize.root(mole_balance, guess)
                        if not(soln.success):
                            print("A solution was NOT obtained:")
                            print(soln.message)
                        CY = soln.x[0]

                        # add +/- 0.01 random "error"
                        random_error = (2*random.random() - 1.0)*0.01
                        CY += random_error

                        # round to 2 decimal places
                        CY = round(CY,2)

                        # append the result to the dataframe
                        df.loc[len(df.index)] = [T,1000*Vdot,CAin,CBin,CYin,CZin,CY]

# display the results
print("\n")
print(df)

# make sure a REB folder exists for the storing results
reb_book_path = '../../../RE_Basics/solutions/ch11_ex2/'
if not os.path.isdir(reb_book_path):
    # create the folder
    os.makedirs(reb_book_path)

# save the results
df.to_csv('example_11_5_2_data.csv',index=False)
df.to_csv(reb_book_path + 'example_11_5_2_data.csv',index=False)
