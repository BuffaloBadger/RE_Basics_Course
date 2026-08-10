"""Calculations for Example 12.4.2 of REB, The Book"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import pandas as pd

# constants available to all functions
# given
dH = -9120 # cal /mol
K0 = 0.132
CpCO = 29.3/4.184 # cal /mol /K
CpH2O = 34.3/4.184 # cal /mol /K
CpCO2 = 41.3/4.184 # cal /mol /K
CpH2 = 29.1/4.184 # cal /mol /K
CpI = 40.5/4.184 # cal /mol /K
nDotCO_0 = 1.0 # mol /h
nDotCO2_0 = 0.359 # mol /h
nDotH2_0 = 4.44 # mol /h
nDotI_0 = 0.18 # mol/h
nDotH2O_0 = 9.32 # mol /h
P = 26.0 # atm
T0 = 445 + 273.15 # K
V1 = 685 # cm^3
k01 = 3.54E-2 # mol /cm^3 /min /atm^2
E1 = 9740 # cal /mol
T4 = 20 + 273.15 # K
mDot4 = 1100 # g /h
T5 = 50 + 273.15 # K
CpEx = 1.0 # cal /g /K
V2 = 3950 # cm^3
k02 = 1.77E-3 # mol /cm^3 /min /atm^2
E2 = 3690 # cal /mol
# known
R = 1.987 # cal /mol /K

# global variable for the reactor number
global g_reactorNo
g_reactorNo = float('nan')

# PFR model function
def pfr_model_variables(dep0, reactorNo):
    # make the reactor number available to the pfr derivatives function
    global g_reactorNo
    g_reactorNo = reactorNo

	# set the initial values
    ind0 = 0.0

	# define the stopping criterion
    fVar = 0
    fVal = V1
    if reactorNo == 2:
        fVal = V2
     
	# solve the design equations
    V, dep, success, message = solve_ivodes(ind0, dep0, fVar, fVal
        , pfr_derivatives, odes_are_stiff=False)

    # check for solver issues
    if not(success):
        print('')
        print(f"PFR model solver issue: {message}")
        print('')
        input('Press return to continue or CTRL-C to exit.')

    # extract the dependent variable pfr_model_variables
    nDotCO = dep[0,:]
    nDotH2O = dep[1,:]
    nDotCO2 = dep[2,:]
    nDotH2 = dep[3,:]
    nDotI = dep[4,:]
    T = dep[5,:]

    # return the pfr_model_variables
    return V, nDotCO, nDotH2O, nDotCO2, nDotH2, nDotI, T

# pfr_derivatives function
def pfr_derivatives(ind, dep):
	# extract the dependent variables
    nDotCO = dep[0]
    nDotH2O = dep[1]
    nDotCO2 = dep[2]
    nDotH2 = dep[3]
    nDotI = dep[4]
    T = dep[5]

	# calculate additional unknowns
    nDotTotal = nDotCO + nDotH2O + nDotCO2 + nDotH2 + nDotI
    PCO = P*nDotCO/nDotTotal
    PH2O = P*nDotH2O/nDotTotal
    PCO2 = P*nDotCO2/nDotTotal
    PH2 = P*nDotH2/nDotTotal
    K = K0*np.exp(-dH/R/T)
    if g_reactorNo == 1:
        k = k01*np.exp(-E1/R/T)
    else:
        k = k02*np.exp(-E2/R/T)
    r = k*(PCO*PH2O - PCO2*PH2/K)
    
	# evaluate the PFR derivatives
    dCOdV = -r
    dH2OdV = -r
    dCO2dV = r
    dH2dV = r
    dIdV = 0
    dTdV = -r*dH/(nDotCO*CpCO + nDotH2O*CpH2O + nDotCO2*CpCO2 \
        + nDotH2*CpH2 + nDotI*CpI)
    
	# return the PFR derivatives
    return dCOdV, dH2OdV, dCO2dV, dH2dV, dIdV, dTdV

# deliverables function
def deliverables():
    # solve the reactor 1 design equations
    reactorNo = 1

    # define the initial values of the dependent variables
    dep0 = np.array([nDotCO_0, nDotH2O_0, nDotCO2_0, nDotH2_0, nDotI_0,
        T0])

    # solve the reactor design equations
    V, nDotCO, nDotH2O, nDotCO2, nDotH2, nDotI, T = pfr_model_variables(dep0
        , reactorNo)

    # extract the outlet temperature and calculate the conversion
    T1 = T[-1]
    fCO_1 = 100*(nDotCO_0 - nDotCO[-1])/nDotCO_0

    # heat exchanger model
    T2 = T1 - mDot4*CpEx*(T5 - T4) \
        /(nDotCO[-1]*CpCO + nDotH2O[-1]*CpH2O \
        + nDotCO2[-1]*CpCO2 + nDotH2[-1]*CpH2 + nDotI[-1]*CpI)

    # solve the reactor 2 design equations
    reactorNo = 2

    # define the initial values of the dependent variables
    dep0 = np.array([nDotCO[-1], nDotH2O[-1], nDotCO2[-1], nDotH2[-1],
        nDotI[-1], T2])

    # solve the reactor design equations
    V, nDotCO, nDotH2O, nDotCO2, nDotH2, nDotI, T = pfr_model_variables(dep0
        , reactorNo)

    # extract the outlet temperature and calculate the conversion
    T3 = T[-1]
    fCO_2 = 100*(nDotCO_0 - nDotCO[-1])/nDotCO_0

    # tabulate, show, and save the results
    results = [['0', 0.0, T0 - 273.15], ['1', fCO_1, T1 - 273.15]
        ,['2', fCO_1, T2 - 273.15], ['3', fCO_2, T3 - 237.15]]
    results_df = pd.DataFrame(results,columns=['stream','conversion'
        ,'temperature'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('results.csv',index=False)
    return

if __name__=="__main__":
    deliverables()
