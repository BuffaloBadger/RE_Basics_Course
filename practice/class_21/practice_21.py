"""Calculations for the Class 21 Practice Assignment from REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import scipy.integrate as sp
from reb_utils import solve_ivodes
from reb_utils import solve_ates
import matplotlib.pyplot as plt

# global constants available to all functions
# given
Tin = 140 + 273.15 # K
P = 1.0 # atm
yAin = 0.3
yBin = 0.7
VdotIn = 450 # cm^3 /min
D = 4 # cm
Tex_in = 140 + 273.15 # K
mDot = 20 # g /min
Cpex = 1.4/1000 # kJ /g /K
Uex = 1.3E-6 # kJ /min /cm^2 /K
tau = 20 # min
k01 = 1.65E4 # mol /min /cm^3 /atm^2
k02 = 3.24E4 # mol /min /cm^3 /atm^2
E1 = 78 # kJ /mol
E2 = 86 # kJ /mol
dH1 = -35.1 # kJ /mol
dH2 = -32.6 # kJ /mol
CpA = 78.3/1000 # kJ /mol /K
CpB = 81.1/1000 # kJ /mol /K
CpD = 75.4/1000 # kJ /mol /K
CpZ = 68.3/1000 # kJ /mol /K
CpU = 76.3/1000 # kJ /mol /K
# known
Re = 8.314/1000 # kJ /mol /K
Rp = 82.057 # cm^3 atm /mol /K
# calculated
nDotAin = yAin*VdotIn*P/Rp/Tin
nDotBin = yBin*VdotIn*P/Rp/Tin
nDotDin = 0
nDotZin = 0
nDotUin = 0
L = tau*VdotIn/((np.pi*D**2/4))

# global variable for the current exchange fluid temperature
g_Tex = float('nan')

# PFR reactor function
def pfr_model_variables(Tex):
    # make Tex available to the derivatives function
    global g_Tex
    g_Tex = Tex

    # define the initial values
    ind_0 = 0
    dep_0 = np.array([nDotAin, nDotBin, nDotDin, nDotZin, nDotUin, Tin])

    # define the stopping criterion
    f_var = 0
    f_val = L

    # solve the design equations
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
        , pfr_derivatives, odes_are_stiff=True)
    
    # check for solver issues
    if not success:
        print('')
        print(f'PFR model issue: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit.')
    
    # return the pfr model variables
    return z, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:]

# PFR derivatives function
def pfr_derivatives(z, dep):
    # extract the dependent variables
    nDotA = dep[0]
    nDotB = dep[1]
    nDotD = dep[2]
    nDotZ = dep[3]
    nDotU = dep[4]
    T = dep[5]

    # calculate the additional unknowns
    nDotTotal = nDotA + nDotB + nDotD + nDotZ + nDotU
    PA = nDotA/nDotTotal*P
    PB = nDotB/nDotTotal*P
    PD = nDotD/nDotTotal*P
    k1 = k01*np.exp(-E1/Re/T)
    k2 = k02*np.exp(-E2/Re/T)
    r1 = k1*PA*PB
    r2 = k2*PD*PB

    # evaluate the derivatives
    dnDotAdz = np.pi*D**2/4*(-r1)
    dnDotBdz = np.pi*D**2/4*(-r1 -r2)
    dnDotDdz = np.pi*D**2/4*(r1 - r2)
    dnDotZdz = np.pi*D**2/4*(r1 + r2)
    dnDotUdz = np.pi*D**2/4*(r2)
    dTdz = (np.pi*D*Uex*(g_Tex-T) - np.pi*D**2/4*(r1*dH1 + r2*dH2))/(nDotA*CpA + nDotB*CpB + nDotD*CpD + nDotZ*CpZ + nDotU*CpU)

    # return the derivatives
    return [dnDotAdz, dnDotBdz, dnDotDdz, dnDotZdz, dnDotUdz, dTdz]

# coupled unknown residual function
def coupled_unknown_residual(guess):
    # solve the pfr model
    z, nDotA, nDotB, nDotD, nDotZ, nDotU, T = pfr_model_variables(guess)

    # calculate Qdot
    Qdot = np.pi*D*Uex*(sp.trapezoid(guess - T, z))

    # evaluate the residual
    epsilon = Qdot + mDot*Cpex*(guess - Tex_in)

    # return the residual
    return epsilon

# deliverables function
def deliverables():
    # solve for the exchange fluid temperature
    guess = Tex_in + 5
    soln, success, message = solve_ates(coupled_unknown_residual, guess)
    
    # check for solver issues
    if not success:
        print('')
        print(f'Coupled unknown issue: {message}')
        print('')
        input('Press return to continue or CTRL-C to exit.')

    # solve the design equations with the solved exchange fluid temperature
    z, nDotA, nDotB, nDotD, nDotZ, nDotU, T = pfr_model_variables(soln[0])

    # calculate the quantities of interest
    conversion = 100*(nDotAin - nDotA[-1])/nDotAin
    selectivity = nDotD[-1]/nDotU[-1]
    T_out = T[-1]
    Tex_out = soln[0]

    # tabulate, show and save the results
    results =[['Conversion', f'{conversion:.2f}','%']
              , ['Selectivity', f'{selectivity:.2f}','mol D per mol U']
              , ['T out', f'{T_out - 273.15:.2f}', '°C']
              , ['Tex out', f'{Tex_out - 273.15:.2f}', '°C']]
    results_df = pd.DataFrame(results, columns=['Item', 'Value', 'Units'])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('practice_21_results.csv', index=False)

    # for discussion, plot the temperature profile
    plt.figure(1)
    plt.plot(z, T - 273.15, label='T')
    plt.axhline(y=soln[0] - 273.15, color='r', linestyle='--', label='Tex')
    plt.xlabel('Axial Position, z (cm)')
    plt.xlim(left=0)
    plt.ylabel('Temperature (°C)')
    plt.legend()
    plt.savefig('practice_21_temperature_profile.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
