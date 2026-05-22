"""Calculations for the Class 19 Learning Activity from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes
from reb_utils import solve_ates

# set resolution for graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
dH_1 = -184000.0 # J/mol
k0_1 = 4.866E14 # L/mol/min
E_1 = 74100.0 # J/mol
P = 1.0 # atm
V0 = 4.0 # L
CA0 = 10.0 # mol /L
T0 = 20 + 273.15 # K
Tex_0 = 20 + 273.15 # K
Vex = 500.0/1.0E3 # L
Aex = 550.0/1.0E4 # m^2
U = 8480.0 # J/m^2/min/K
Vdot_ex = 1.0 # L/min
Tex0 = 20 + 273.15 # K
Cp = 4.184 # J/g/K
rho = 1.0E3 # g /L
CBin = 10.0 # mol/L
Tex_in = 20 + 273.15 # K
Tin = 20 + 273.15 # K
Tmax = 80.0 + 273.15 # K
# known
R = 8.3136 # J/mol/K
Rpv = 8.206E-2 # L-atm/mol/K
# calculated
nA0 = CA0*V0
nAf = 0.01*nA0

# global variable for the current feed rate
g_Vdot = float('nan')

# SBSTR reactor function
def sbstr_model_variables(Vdot):
    # make Vdot available to the derivatives function
    global g_Vdot
    g_Vdot = Vdot

    # define the initial values
    ind_0 = 0
    dep_0 = np.array([nA0, 0, 0, T0, Tex0, V0])

    # define the stopping criterion
    f_var = 1
    f_val = nAf

    # solve the SBSTR design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            ,sbstr_derivatives, odes_are_stiff=True)
    
    # check for solver issues
    if not success:
        print('')
        print(f'SBSTR model function issue: {message}')
        print('')
        input('Press return to continue.')
    
    # return the sbstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:]

# SBSTR derivatives function
def sbstr_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nZ = dep[2]
    T = dep[3]
    Tex = dep[4]
    V = dep[5]

    # calculate the additional unknowns
    nBin = g_Vdot*CBin
    k1 = k0_1*np.exp(-E_1/R/T)
    CA = nA/V
    CB = nB/V
    r1 = k1*CA*CB
    Qdot = U*Aex*(Tex - T)
    mDot_ex = Vdot_ex*rho

    # evaluate the derivatives
    dnAdt = -r1*V
    dnBdt = nBin - r1*V
    dnZdt = 2*r1*V
    dTdt = (Qdot - g_Vdot*Cp*rho*(T - Tin) - r1*V*dH_1 + P*g_Vdot*R/Rpv)/(V*Cp*rho)
    dTexdt = -(Qdot + mDot_ex*Cp*(Tex - Tex_in))/(rho*Vex*Cp)
    dVdt = g_Vdot

    # return the derivatives
    return np.array([dnAdt, dnBdt, dnZdt, dTdt, dTexdt, dVdt])

# coupled unknown residual function
def coupled_unknown_residual(guess):
    # solve the sbstr design equations using the guess
    t, nA, nB, nZ, T, Tex, V = sbstr_model_variables(guess[0])

    # evaluate and return the residual
    epsilon = max(T) - Tmax
    return epsilon

# deliverables function
def deliverables():
    # guess the optimum flow rate
    #Vdot_guess = 1.0.  failed to converge
    Vdot_guess = 0.01

    # solve the implicit equation for the optimum flow rate
    soln, success, message = solve_ates(coupled_unknown_residual, Vdot_guess)
    Vdot_opt = soln[0]

    # check for solver issues
    if not success:
        print('')
        print(f'Issue solving the implicit ATE: {message}')
        print('')
        input('Press return to continue.')
    
    # solve the sbstr design equations using the optimum Vdot
    t, nA, nB, nZ, T, Tex, V = sbstr_model_variables(Vdot_opt)

    # find the maximum T and concentration of B
    maxT = max(T) - 273.15
    max_CB = max(nB/V)

    # tabulate, show, and save the results
    data = [["Optimum Feed Rate", Vdot_opt, "L/min"]
            ,["Maximum Temperature", maxT, "°C"]
            ,["Maximum Concentration of B", max_CB, "M"]]
    results_df = pd.DataFrame(data, columns = ["Item", "Value", "Units"])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_19_results.csv', index=False)

    # for discussion, plot fA and T vs t
    plt.figure(1)
    plt.plot(t,100*(nA0 - nA)/nA0)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Conversion of A (%)')
    plt.ylim(bottom=0, top=100)
    plt.savefig('activity_19_fA_vs_t.pdf')

    plt.figure(2)
    plt.plot(t, T-273.15, label='Reacting Fluid')
    plt.plot(t,Tex-273.15, label='Cooling Water')
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Temperature (°C)')
    plt.legend()
    plt.savefig('activity_19_T_vs_t.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    