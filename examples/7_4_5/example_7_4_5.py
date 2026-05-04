"""Calculations for Example 7.4.5 from REB, The Book"""

# import libraries
import numpy as np
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes

# set resolution of graphs
plt.rc('savefig', dpi=300)

# global constants available to all functions
# given
k01 = 1.6E18 # gal /lbmol /h
E1 = 46000. # BTU /lbmol
k02 = 4.5E18 # gal /lbmol /h
E2 = 48000. # BTU /lbmol
Cp = 65.0 * 0.1337 # BTU /gal /degR
dH1 = 45000. # BTU /lbmol
dH2 = 39500. # BTU /lbmol
V = 50 # gal
CA0 = 0.014 # lbmolA /gal
CB0 = 0.020 # lbmolB /gal
T0 = 70 + 459.67 # degR
U = 60 # BTU /ft^2 /degR /h
A = 13 # ft^2
Tex = 220 + 459.67 # degR
dHvap = 966 # BTU/lb
# known
R = 1.987 # BTU/lbmol
# calculated
nA0 = CA0*V
nB0 = CB0*V

# BSTR model function
def bstr_model_variables(fB):
    # set the initial values
    ind_0 = 0
    dep_0 = np.array([nA0, nB0, 0, 0, 0, T0])

    # set the stopping criterion
    f_var = 2
    f_val = nB0*(1-fB)

     # solve the design equations
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                ,bstr_derivatives, odes_are_stiff=False)
    
    # check for solver issues
    if not success:
        print('')
        print(f"BSTR model function issue: {message}")
        print('')
        input("Press return to continue.")
    
    # return the bstr model variables
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:]

# BSTR derivatives function
def bstr_derivatives(ind, dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]
    nD = dep[2]
    nZ = dep[3]
    nU = dep[4]
    T = dep[5]

    # calculate the additional unknowns
    k1 = k01*np.exp(-E1/R/T)
    k2 = k02*np.exp(-E2/R/T)
    CA = nA/V
    CB = nB/V
    CD = nD/V
    r1 = k1*CA*CB
    r2 = k2*CD*CB
    Q = U*A*(Tex-T)

    # evaluate the derivatives
    dnAdt = -r1*V
    dnBdt = -(r1 + r2)*V
    dnDdt = (r1 - r2)*V
    dnZdt = (r1 + r2)*V
    dnUdt = r2*V
    dTdt = (Q - (r1*dH1 + r2*dH2)*V)/V/Cp

    # return the derivatives
    return np.array([dnAdt, dnBdt, dnDdt, dnZdt, dnUdt, dTdt])

# deliverables function
def deliverables():
    # solve the BSTR design equations
    t, nA, nB, nD, nZ, nU, T = bstr_model_variables(0.99)

    # calculate the quantities to be plotted
    fB = 100*(nB0 - nB)/nB0
    Yield = nD/nA0
    Qdot = U*A*(Tex - T)
    m_H2O = Qdot/dHvap

    # generate, show, and save the graphs
    plt.figure(1)
    plt.plot(fB,Yield)
    plt.ylabel('Yield (mol D per mol A fed)')
    plt.xlabel('Conversion of B (%)')
    plt.xlim(left=0, right=100)
    plt.ylim(bottom=0)
    plt.savefig('example_7_4_5_yield_vs_conversion.png')
    plt.savefig('example_7_4_5_yield_vs_conversion.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex5/example_7_4_5_yield_vs_conversion.png')
    plt.show(block=False)

    plt.figure(2)
    plt.plot(fB,T - 459.67)
    plt.ylabel('Temperature (°F)')
    plt.xlabel('Conversion of B (%)')
    plt.xlim(left=0, right=100)
    plt.savefig('example_7_4_5_T_vs_conversion.png')
    plt.savefig('example_7_4_5_T_vs_conversion.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex5/example_7_4_5_T_vs_conversion.png')
    plt.show(block=False)

    plt.figure(3)
    plt.plot(fB,m_H2O)
    plt.ylabel('Water Flow Rate (lb$_m$ h$^{-1}$)')
    plt.xlabel('Conversion of B (%)')
    plt.xlim(left=0, right=100)
    plt.savefig('example_7_4_5_H2O_vs_conversion.png')
    plt.savefig('example_7_4_5_H2O_vs_conversion.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex5/example_7_4_5_H2O_vs_conversion.png')
    plt.show()

    return

# execution command
if __name__ == '__main__':
    deliverables()
