"""Calculations for the Class 15 Learning Activity from REB, The Book"""

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
k01 = 2.59e9 # min-1
E1 = 16500 # cal/mol
dH1 = -22200 # cal/mol
Cp = 440 # cal/L/K
V = 4 # L
CA_charge = 2 # mol/l
T0 = 60 + 273.15 #K
fA0 = 0.18
Tex0 = 60 + 273.15 #K
Tex_in = 60 + 273.15 # K
mDotEx = 1.5 # kg/min
rhoEx = 1 # kg/L
Cpex = 1000 # cal/kg/K
CAf = 0.2 # mol/L
Vex = 0.5 #L
U = 1.13e4/60 # cal/ft2/min/K
tf = 20 # min
# known
R = 1.987 # cal/mol/K
# calculated
nA0 = CA_charge*V*(1-fA0)
nZ0 = CA_charge*V*fA0
nAf = CAf*V

# Allocate storage to make Aex globally available
g_Aex = float('NaN')

# BSTR reactor function
def bstr_model_variables(Aex):
    # make Aex available to the bstr residuals function
    global g_Aex
    g_Aex = Aex

    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nA0, nZ0, T0, Tex0])

    # set the stopping criterion
    f_var = 0
    f_val = tf

    # solve the design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
            ,bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model function issue: {message}")
        print('')
        input('Press return to continue.')

    # return the bstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:]

# BSTR derivatives function
def bstr_derivatives(ind,dep):
    # extract the dependent variables
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]
    Tex = dep[3]

    # calculate the additional unknowns
    k1 = k01*np.exp(-E1/(R*T))
    CA = nA/V
    r1 = k1*CA
    Qdot = U*g_Aex*(Tex - T)

    # evaluate the derivatives
    dnAdt = -r1*V
    dnZdt = r1*V
    dTdt = (Qdot -V*r1*dH1)/(V*Cp)
    dTexdt = (mDotEx*Cpex*(Tex_in - Tex) - Qdot)/(rhoEx*Vex*Cpex)

    # return the derivatives
    return dnAdt, dnZdt, dTdt, dTexdt

# coupled unknown residual function
def coupled_unknown_residual(AexGuess):

    # solve the BSTR design equations
    t, nA, nZ, T, Tex = bstr_model_variables(AexGuess)

    # evaluate the residual
    CA_fromGuess = nA[-1]/V
    epsilon = CAf - CA_fromGuess

    # return the residual
    return epsilon

# deliverables function
def deliverables():
    # guess Aex
    Aex = 1 # ft^2

    # calculate Aex
    soln, success, message = solve_ates(coupled_unknown_residual,Aex)
    Aex = soln[0]

    # check for solver issues
    if not success:
        print('')
        print(f"Issue solving for the coupled unknown: {message}")
        print('')
        input('Press return to continue.')

    # solve the BSTR design equations
    t, nA, nZ, T, Tex = bstr_model_variables(Aex)

    # tabulate, show, and save the deliverables
    data = [["Area",f"{Aex:.2f}","square feet"]
            ,["Reacting Fluid Temperature",f"{T[-1]-273.15:.1f}","°C"]]
    results_df = pd.DataFrame(data,columns=["Item", "Value", "Units"])
    print('')
    print(results_df)
    print('')
    results_df.to_csv('activity_15_results.csv',index=False)

    # for discussion, plot T vs t
    plt.figure(1)
    plt.plot(t, T-273.15)
    plt.xlabel('Time (min)')
    plt.xlim(left=0)
    plt.ylabel('Reacting Fluid Temperature (°C)')
    plt.savefig('activity_15_T_vs_t.png')
    plt.savefig('activity_15_T_vs_t.pdf')
    plt.show()

# execution command
if __name__ == '__main__':
    deliverables()
    