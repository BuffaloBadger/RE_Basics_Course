"""Calculations for modification of the 2021 Assessment #5 for use as the class 21 practice assignment."""

# import libraries
import numpy as np
import pandas as pd
import scipy.integrate as sp
import matplotlib.pyplot as plt
from reb_utils import solve_ates
from reb_utils import solve_ivodes

# given or known
Re = 8.314/1000 # kJ /mol /K
Rp = 82.057 # cm^3 atm /mol /K
#Tin = 300 + 273.15 # K
P = 1.0 # atm
yAin = 0.3
yBin = 0.7
VdotIn = 450 # cm^3 /min
D = 4 # cm
#Te = 160 + 273.15 # K
#U = 0.0013 # kJ /min /cm^2 /K
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
nDotDin = 0
nDotZin = 0
nDotUin = 0

Tin = 140 + 273.15
Te = 145 + 273.15 # K
tau = 20 # min
L = tau*VdotIn/((np.pi*D**2/4))
U = 1.3e-6
mDot = 20 # g /min
Cpex = 1.4/1000

# Calculated constants
nDotAin = yAin*VdotIn*P/Rp/Tin
nDotBin = yBin*VdotIn*P/Rp/Tin

# PFR model function
def pfr_model_variables():
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
    dTdz = (np.pi*D*U*(Te-T) - np.pi*D**2/4*(r1*dH1 + r2*dH2))/(nDotA*CpA + nDotB*CpB + nDotD*CpD + nDotZ*CpZ + nDotU*CpU)

    # return the derivatives
    return np.array([dnDotAdz, dnDotBdz, dnDotDdz, dnDotZdz, dnDotUdz, dTdz])

def deliverables():
    # solve the design equations
    z, nDotA, nDotB, nDotD, nDotZ, nDotU, T = pfr_model_variables()

    # calculate Qdot
    Qdot = np.pi*D*U*(sp.trapezoid(Te - T, z))
    TexIn = (Qdot + mDot*Cpex*Te)/(mDot*Cpex)
    yld = 100*nDotD[-1]/nDotAin
    sel = nDotD[-1]/nDotU[-1]

    # display the results
    print('')
    print(f'Qdot = {Qdot:.2f} kJ/min')
    print(f'TexIn = {TexIn - 273.15:.2f} K')
    print(f'Yield = {yld:.2f} %')
    print(f'Selectivity = {sel:.2f} mol D per mol U')
    print('')

    plt.figure(1)
    plt.plot(z, T - 273.15, label='T')
    plt.xlabel('Position (cm)')
    plt.ylabel('Temperature (°C)')
    plt.title('PFR Temperature Profile')
    plt.legend()
    plt.show(block=False)

    plt.figure(2)
    plt.plot(z, nDotA, label='nDotA')
    plt.plot(z, nDotB, label='nDotB')
    plt.plot(z, nDotD, label='nDotD')
    plt.plot(z, nDotZ, label='nDotZ')
    plt.plot(z, nDotU, label='nDotU')
    plt.xlabel('Position (cm)')
    plt.ylabel('Molar Flow Rate (mol/min)')
    plt.title('PFR Molar Flow Rate Profile')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    deliverables()