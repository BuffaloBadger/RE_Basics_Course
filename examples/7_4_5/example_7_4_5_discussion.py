"""Calculations for discussion of Example 7.4.5 from REB, The Book"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from reb_utils import solve_ivodes
import example_7_4_5

# set resolution of graphs
plt.rc('savefig', dpi=300)

# deliverables function
def deliverables():
    # solve the BSTR design equations
    t, nA, nB, nD, nZ, nU, T = example_7_4_5.bstr_model_variables(0.99)

    # calculate the quantities to be plotted
    fB = 100*(example_7_4_5.nB0 - nB)/example_7_4_5.nB0
    Yield = nD/example_7_4_5.nA0
    Qdot = example_7_4_5.U*example_7_4_5.A*(example_7_4_5.Tex - T)
    m_H2O = Qdot/example_7_4_5.dHvap

    plt.figure(1)
    plt.plot(t,fB)
    plt.ylabel('Conversion of B (%)')
    plt.xlabel('Time (h)')
    plt.xlim(left=0)
    plt.ylim(bottom=0, top=100)
    plt.savefig('example_7_4_5_fB_vs_t.png')
    plt.savefig('example_7_4_5_fB_vs_t.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex5/example_7_4_5_fB_vs_t.png')
    plt.show(block=False)

    # generate, show, and save the graphs
    plt.figure(2)
    plt.plot(t,Yield)
    plt.ylabel('Yield (mol D per mol A fed)')
    plt.xlabel('Time (h)')
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.savefig('example_7_4_5_yield_vs_t.png')
    plt.savefig('example_7_4_5_yield_vs_t.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex5/example_7_4_5_yield_vs_t.png')
    plt.show(block=False)

    plt.figure(3)
    plt.plot(t,T - 459.67)
    plt.ylabel('Temperature (°F)')
    plt.xlabel('Time (h)')
    plt.xlim(left=0)
    plt.savefig('example_7_4_5_T_vs_t.png')
    plt.savefig('example_7_4_5_T_vs_t.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex5/example_7_4_5_T_vs_t.png')
    plt.show(block=False)

    plt.figure(4)
    plt.plot(t,m_H2O)
    plt.ylabel('Water Flow Rate (lb$_m$ h$^{-1}$)')
    plt.xlabel('Time (h)')
    plt.xlim(left=0)
    plt.savefig('example_7_4_5_H2O_vs_t.png')
    plt.savefig('example_7_4_5_H2O_vs_t.pdf')
    plt.savefig('../../../RE_Basics/solutions/ch7_ex5/example_7_4_5_H2O_vs_t.png')
    plt.show()

    return

# execution command
if __name__ == '__main__':
    deliverables()
