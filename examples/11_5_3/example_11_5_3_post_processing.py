"""Post processing for Example 11.5.3 from REB, The Book"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from great_tables import GT, md, html

# make sure a png folder exists
if not os.path.isdir('png'):
    # create the folder
    os.makedirs('./png')

# table with the first six rows of the experimental data
table_df = pd.read_csv('example_11_5_3_data.csv')
table_df = table_df.rename(columns={'T (F)': 'T', 'Vdot (gal/min)': 'Vdot', 'CA_in (lbmol/gal)': 'CA_in', 'CY_in (lbmol/gal)': 'CY_in', 'CZ_in (lbmol/gal)': 'CZ_in', 'CA_out (lbmol/gal)': 'CA_out'})
table_df = table_df.head(6)
data_tbl = (GT(table_df)
    .cols_label(T=html('T <br> (°F)'),
                Vdot=html("""<span style="position:relative;">    V<span style="position:absolute; left:0.1em; top:-0.3em;">˙</span></span><br>(gal min<sup>-1</sup>)"""),
                CA_in=html('C<sub>A,in</sub> <br> (lbmol gal<sup>-1</sup>)'),
                CY_in=html('C<sub>Y,in</sub> <br> (lbmol gal<sup>-1</sup>)'),
                CZ_in=html('C<sub>Z,in</sub> <br> (lbmol gal<sup>-1</sup>)'),
                CA_out=html('C<sub>Y,out</sub> <br> (lbmol gal<sup>-1</sup>)')
    )
    .fmt_number(columns=[0,1], decimals=0)
    .cols_align(
        align='center',
        columns=[0,1,2,3,4,5]
    )
)
#data_tbl.show()
data_tbl.gtsave('png/data_table.png', zoom=4.0)