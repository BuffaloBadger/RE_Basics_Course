"""Calculations for discussion of Example 6.6.4 of REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ates

# set the dpi for figures
plt.rc("savefig", dpi=300)

# constants available to all functions
# given
V = 500 # cm^3
Vdot_in = 1.0 # cm^3 /s
CA_in = 0.015 # mol /cm^3
CB_in = 0.015 # mol /cm^3
T_in = 50 # °C
T_s_hi = 264.75 # °C
T_s_lo = 50.06 # °C
T_u = 137.96 # °C
Cp = 0.35*4.184 # J /g /K
rho = 0.93 # g /cm^3
dH = -20000 # J /mol
k0 = 3.24E12 # cm^3 /mol /s
E = 105000 # J /mol
# known
R = 8.314 # J /mol /K
# calculated
nA_in = Vdot_in*CA_in
nB_in = Vdot_in*CB_in

# make T available to all functions
T = float('NaN')

# reactor model function
def reactor_model_variables(init_guess):
     
	# solve the ATEs
    soln, success, message = solve_ates(residuals,init_guess)

    # check that the solution is converged
    if not(success):
        print(f"solve_ates was unsuccessful: {message}")

    # return the solution
    return soln

# residuals function
def residuals(guess):
    # extract the individual guesses
    nA = guess[0]
    nB = guess[1]
    nY = guess[2]
    nZ = guess[3]
    T_in = guess[4]

    # calculate the rate
    k = k0*np.exp(-E/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r = k*CA*CB

    # evaluate the residuals
    residual_1 = nA_in - nA - V*r
    residual_2 = nB_in - nB - V*r
    residual_3 = -nY + V*r
    residual_4 = -nZ + V*r
    residual_5 = -Vdot_in*rho*Cp*(T - T_in) - V*r*dH

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4
                     , residual_5])

# deliverables function
def deliverables():

    # enable this function to set T
    global T
    
    # set a range of outlet temperatures
    T_range = np.linspace(-25.0, 325.0, 350) + 273.15

	# set an initial guess
    init_guess = np.array([0.99*nA_in, 0.99*nB_in, 0.01*nA_in, 0.01*nA_in
            ,T_range[0] - 1.0])

    # allocate storage for the conversions and inlet temperatures
    f_range = np.zeros(350)
    Tin_range = np.zeros(350)

    # calculate the conversion and inlet temperature for each outlet temperature
    for i in range(0,350):

        # set the outlet temperature
        T = T_range[i]

        # solve the design equations
        solution = reactor_model_variables(init_guess)

        # save the conversion and inlet temperature
        nA = solution[0]
        f_range[i] = 100*(nA_in - nA)/nA_in
        Tin_range[i] = solution[4]

        # use the result as the initial guess for the next T
        init_guess = solution

    # convert to °C
    Tin_range = Tin_range - 273.15
    T_range = T_range - 273.15

    # find the extinction and ignition points
    iLast = 0
    need_ignition = True
    for i in range(1,350):
        if need_ignition:
            if Tin_range[i] < Tin_range[iLast]:
                i_ignition = iLast
                need_ignition = False
            iLast = i
        else:
            if Tin_range[i] > Tin_range[iLast]:
                i_extinction = iLast
                break
            iLast = i
   
    # tabulate, show, and save the ignition and extinction points
    data = [['Extinction Point Feed T',f'{Tin_range[i_extinction]}','°C']
            ,['Extinction Point T',f'{T_range[i_extinction]}','°C']
            ,['Extinction Point Conversion',f'{f_range[i_extinction]}','%']
            ,['Ignition Point Feed T',f'{Tin_range[i_ignition]}','°C']
            ,['Ignition Point T',f'{T_range[i_ignition]}','°C']
            ,['Ignition Point Conversion',f'{f_range[i_ignition]}','%']]
    results_df = pd.DataFrame(data, columns=['item' ,'value','units'])
    print(results_df)
    results_df.to_csv('example_6_6_4_ign_ext_points.csv',index=False)

    # multiplicity plot for section 6.4
    plt.figure(2) 
    plt.plot(Tin_range[0:i_ignition], T_range[0:i_ignition],'k-')
    plt.plot(Tin_range[i_ignition+1:i_extinction]
             , T_range[i_ignition+1:i_extinction],'k--')
    plt.plot(Tin_range[i_extinction+1:-1], T_range[i_extinction+1:-1],'k-')
    plt.axvline(T_in, c='b')
    plt.axvline(Tin_range[i_extinction],c='r',ls='-')
    plt.axvline(Tin_range[i_ignition],c='r',ls='-')
    plt.plot(T_in, T_s_hi, 'ko')
    plt.text(T_in + 3,T_s_hi - 20, "$S_{high}$")
    plt.plot(T_in,T_u, 'ko')
    plt.text(T_in + 3,T_u + 10, "$U_{mid}$")
    plt.plot(T_in,T_s_lo, 'ko')
    plt.text(T_in + 3,T_s_lo - 20, "$S_{low}$")
    plt.xlabel("Inlet CSTR Temperature")
    plt.ylabel("Outlet CSTR Temperature")
    plt.xticks([])
    plt.yticks([])
    plt.savefig('CSTR_multiplicity.png')
    plt.savefig('../../../RE_Basics/figures/CSTR_multiplicity.png')
    plt.show()
    return

# calculate the deliverables
if __name__=="__main__":
    deliverables()
