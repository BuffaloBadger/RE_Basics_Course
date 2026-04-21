"""Calculations for the Class 10 Practice Assignment for REB, The Course"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ates

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants
# given
V = 50 # gal
CA_in = 0.014 # lbmolA /gal
CB_in = 0.020 # lbmolB /gal
T_in = 70 + 459.67 # °R
Cp = 65.0 * 0.1337 # BTU /gal /°R
Tex = 220 + 459.67 # °R
U = 60 # BTU /ft^2 /°R /h
A = 13 # ft^2
dH1 = 45000. # BTU /lbmol
dH2 = 39500. # BTU /lbmol
k01 = 1.6E18 # gal /lbmol /h
E1 = 46000. # BTU /lbmol
k02 = 4.5E18 # gal /lbmol /h
E2 = 54000. # BTU /lbmol
# known
R = 1.987 # BTU/lbmol
# calculated

# allocate global storage for the current value of Vdot_in
g_Vdot_in = float('nan')

# CSTR reactor model function
def cstr_model_variables(Vdot_in, guess):
    # make Vdot_in available to the residuals function
    global g_Vdot_in
    g_Vdot_in = Vdot_in

    # solve the cstr design equations
    soln, success, message = solve_ates(cstr_residuals, guess)

    # check for solver issues
    if not success:
        print("  CSTR model function error: " + message)
    
    # return the solution
    return soln

# CSTR residuals function
def cstr_residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nB = guess[1]
    nD = guess[2]
    nU = guess[3]
    nZ = guess[4]
    T = guess[5]

    # calculate the additional unknowns
    nA_in = CA_in * g_Vdot_in
    nB_in = CB_in * g_Vdot_in
    k1 = k01 * np.exp(-E1 / (R * T))
    k2 = k02 * np.exp(-E2 / (R * T))
    Vdot = g_Vdot_in
    CA = nA / Vdot
    CB = nB / Vdot
    CD = nD / Vdot
    r1 = k1 * CA * CB
    r2 = k2 * CD * CB
    Qdot = U * A * (Tex - T)

    # evaluate the residuals
    epsilon_1 = nA_in - nA - r1 * V
    epsilon_2 = nB_in - nB - r1 * V - r2 * V
    epsilon_3 = -nD + r1 * V - r2 * V
    epsilon_4 = -nU + r2 * V
    epsilon_5 = -nZ + r1*V + r2*V
    epsilon_6 = g_Vdot_in*Cp*(T - T_in) - Qdot + r1*V*dH1 + r2*V*dH2

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4, epsilon_5, epsilon_6])

# deliverables function
def deliverables():
    # choose a range of values for Vdot_in
    VFR_range = np.linspace(35, 50)
    # previous ranges
    #VFR_range = np.linspace(10, 100)
    #VFR_range = np.linspace(20, 40)

    # allocate storage for the corresponding values of nD
    nD_range = np.ones_like(VFR_range) * float('nan')

    # define an initial guess for the first Vdot_in in the range
    guess = np.array([CA_in*VFR_range[0], CB_in*VFR_range[0], 0.0, 0.0, 0.0, T_in + 5])

    # loop through the range of Vdot_in values
    for i, Vdot_in in enumerate(VFR_range):
        # solve the cstr design equations
        soln = cstr_model_variables(Vdot_in, guess)

        # save nD
        nD_range[i] = soln[2]

        # save the solution to use as the next guess
        guess = soln
    
    # plot the outlet molar flow of D versus the inlet volumetric flow rate
    plt.figure(1)
    plt.plot(VFR_range, nD_range)
    plt.ylabel('Outlet Molar Flow Rate of D (lbmol h$^{-1}$)')
    plt.xlabel("Inlet Volumetric Flow Rate (gal h$^{-1}$)")
    plt.savefig('practice_10_nD_vs_VdotIn.png')
    plt.savefig('practice_10_nD_vs_VdotIn.pdf')
    plt.show()

    # find the optimum Vdot_in
    i_max = np.argmax(nD_range)
    Vdot_opt = VFR_range[i_max]

    # solve the cstr design equations at the optimum Vdot_in
    soln = cstr_model_variables(Vdot_opt, guess)

    # calculate the conversion
    nA_in = CA_in*Vdot_opt
    fA = 100*(nA_in - soln[0])/nA_in

    # tabulate, show, and save the results
    data = [['Optimum Vdot_in',Vdot_opt,'gal/h'],
            ['Outlet T',soln[5] - 459.67,'°F'],
            ['Conversion', fA, '%']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])
    print('')
    print(results_df)
    results_df.to_csv('practice_10_results.csv',index = False)

# execution command
if __name__ == "__main__":
    deliverables()
